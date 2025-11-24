#Script used for the transcriptome analysis of Osedax japonicus 
#Results published in "Predetermined sex revealed by a female transient gut in non-feeding larvae of Osedax (Siboglinidae, Annelida)"
#Authors: Alice Rouan, Norio Miyamoto, Katrine Worsaae
#10.1186/s13227-025-00251-9
#https://evodevojournal.biomedcentral.com/articles/10.1186/s13227-025-00251-9
#Gut cluster and transcriptome analysis 
#Alice Rouan's script in K. Worsaae's lab 2025

#Replace path by your own path
library(datapasta)
library(readr)
library(readxl)
library(tidyverse)
#Load data
targets <- read.csv2("path/sample_design.csv", sep=",")
log2.cpm.filtered.norm.df <- read.csv2( "path/Ojap_log2.cpm.filtered.df.csv", row.names="geneID" )
#remove the first column (only row.names)
log2.cpm.filtered.norm.df <-log2.cpm.filtered.norm.df[,-1]
colnames(log2.cpm.filtered.norm.df) <- paste(targets$sample)

#annotation
annotation <- read.csv2("path/Supplementary_material_Table S5_RNAseq_annotation.csv")

#DGE list
myDGEList.filtered.norm <- read.table( "path/myDGEList.filtered.norm.tsv")
#rename the columns
colnames(myDGEList.filtered.norm) <- paste(targets$sample)

##Volcano plot and FC
# Load packages -----
library(tidyverse) # you know it well by now!
library(limma) # venerable package for differential gene expression using linear modeling
library(edgeR)
library(gt) 
library(DT) 
library(plotly) 
library(ggplot2)
library(RColorBrewer)
library(tidyverse) # provides access to Hadley Wickham's collection of R packages for data science, which we will use throughout the course
library(tximport) # package for getting Kallisto results into R
library(datapasta)
library(readr)
library(readxl)

#Group your sample depending on their stage and sex
group <- factor(targets$group_sex)
design <- model.matrix(~0 + group)
colnames(design) <- levels(group)
# fit a linear model to your data
v.DEGList.filtered.norm <- voom(myDGEList.filtered.norm, design, plot = TRUE)
fit <- lmFit(v.DEGList.filtered.norm, design)

#Make a contrast matrix comparing the samples with a gut 4D, 5D and juvenile females vs same stages without a gut Males
contrast.matrixGut <- makeContrasts(GutFvsNoGutM = (D4Larva_F+D5Larva_F+juvenile_F) - (D4Larva_M+D5Larva_M+juvenile_M),levels=design)
fitsGut <- contrasts.fit(fit, contrast.matrixGut)
ebFitGut <- eBayes(fitsGut)

#Now use a for loop 
#More useful when you have multiple comparison
#Chaneg the "path" in the loop
ebFit <- ebFitGut
comparison <- colnames(contrast.matrixGut)

for (i in 1:ncol(ebFit)){
  
  if (colnames(ebFit[,i]$contrasts) == comparison[i])  { #make sure you are using the good comparison
    
    #number=12000 gives a very high cutoff og LogFC 
    myTopHits1 <- topTable(ebFit[,i], adjust ="BH", coef=1, number=140000, sort.by="logFC")

    
    # convert to a tibble
    myTopHits1.df <- myTopHits1 %>%
      as_tibble(rownames = "geneID")

    
    #Add a column to annotate the down,up, non significant or low FC transcript 
    #(up,down,ns,lowFC)
    #for easier colouring when plotting
    #Or even for later on analysis
    
    myTopHits1.df$Cutoff <- NA
    
    myTopHits1.df[which(myTopHits1.df$logFC>0),"Cutoff"] <- "Up"
    myTopHits1.df[which(myTopHits1.df$logFC<0),"Cutoff"] <- "Down"
    myTopHits1.df[which(myTopHits1.df$adj.P.Val >= 0.01),"Cutoff"] <- "Ns"
    
    if(length(which(myTopHits1.df$logFC <= 2 & myTopHits1.df$logFC >= -2 )) == 0){
      
      
      label <- as.data.frame(do.call(rbind,str_split(myTopHits1.df$geneID,"[_]")))
      label[label$V1 %in% "TRINITY", 'V1'] <- NA
      myTopHits1.df$geneLab <-label$V1
      
      #FC cutoff = 2 
      vplot <- ggplot(myTopHits1.df) +
        aes(y=-log10(adj.P.Val), x=logFC, text = paste("Symbol:", geneID)) +
        geom_point(data=myTopHits1.df, aes(colour = factor(Cutoff),fill=factor(Cutoff), alpha=0.2)) +
        #levels Cutoff Down LowFC Ns Up
        scale_colour_manual(values=c("#fde0dd","black","#b8e186"))+
        scale_fill_manual(values=c("#fde0dd","black","#b8e186"))+
        geom_hline(yintercept = -log10(0.01), linetype="longdash", colour="grey", size=0.5) +
        geom_vline(xintercept = 2, linetype="longdash", colour="black", size=0.5) +
        geom_vline(xintercept = -2, linetype="longdash", colour="black", size=0.5) +
             geom_text(aes(label = geneLab, y=-log10(adj.P.Val)+0.5, x=logFC), size =1.8, check_overlap = T)+
        labs(title="Volcano plot",
             subtitle = paste(comparison[i]),
             #caption=paste0("produced on ", Sys.time())) +
             caption=paste0("FC cutoff=2, adj.P.-value cutoff=0.01 ")) +
        theme(text = element_text(size=40))+
        theme_bw()
      
      ggsave(plot = vplot, filename= paste("path/VolcanoPlot_",comparison[i],"_genelevel.pdf"), dpi = 600, width = 15, height = 10)
      ggsave(plot = vplot, filename= paste("path/VolcanoPlot_",comparison[i],"_genelevel.jpg"), dpi = 600, width = 15, height = 10)
      
      write.csv2(myTopHits1.df,file=paste("path/DEG_TopTable_",comparison[i],".csv"))
      
    } else {
      myTopHits1.df[which(myTopHits1.df$logFC <= 2 & myTopHits1.df$logFC >= -2 ),"Cutoff"] <- "LowFC"
      
      
      
      label <- as.data.frame(do.call(rbind,str_split(myTopHits1.df$geneID,"[_]")))
      label[label$V1 %in% "TRINITY", 'V1'] <- NA
      myTopHits1.df$geneLab <-label$V1
      
      #FC cutoff = 2 
      vplot <- ggplot(myTopHits1.df) +
        aes(y=-log10(adj.P.Val), x=logFC, text = paste("Symbol:", geneID)) +
        geom_point(data=myTopHits1.df, aes(colour = factor(Cutoff),fill=factor(Cutoff))) +
        #levels Cutoff Down LowFC Ns Up
        scale_colour_manual(values=c("#fde0dd","black","grey","#b8e186"))+
        scale_fill_manual(values=c("#fde0dd","black","grey","#b8e186"))+
        geom_hline(yintercept = -log10(0.01), linetype="longdash", colour="grey", size=0.5) +
        geom_vline(xintercept = 2, linetype="longdash", colour="black", size=0.5) +
        geom_vline(xintercept = -2, linetype="longdash", colour="black", size=0.5) +
        #annotate("rect", xmin = 1, xmax = 12, ymin = -log10(0.01), ymax = 7.5, alpha=.2, fill="#BE684D") +
        #annotate("rect", xmin = -1, xmax = -12, ymin = -log10(0.01), ymax = 7.5, alpha=.2, fill="#2C467A") +
        geom_text(aes(label = geneLab, y=-log10(adj.P.Val)+0.5, x=logFC), size =1.8, check_overlap = T)+
        labs(title="Volcano plot",
             subtitle = paste(comparison[i]),
             #caption=paste0("produced on ", Sys.time())) +
             caption=paste0("FC cutoff=2, adj.P.-value cutoff=0.01 ")) +
        theme(text = element_text(size=60))+
        theme_bw()
      
      
      ggsave(plot = vplot, filename= paste("path/VolcanoPlot_",comparison[i],"_genelevel.pdf"),dpi = 600, width = 15, height = 10)
      ggsave(plot = vplot, filename= paste("path/VolcanoPlot_",comparison[i],"_genelevel.jpg"),dpi = 600, width = 15, height = 10)
      
      write.csv2(myTopHits1.df,file=paste("path/DEG_TopTable_",comparison[i],".csv"))
      
    }
    
    
  } 
}

#Get the tables generated above
DEG_GutFvsNoGutM <- read.csv2("path/DEG_TopTable_ GutFvsNoGutM .csv")
DEG_GutFvsNoGutM<- DEG_GutFvsNoGutM[,-1]

#Now let's split into Down regulated and Up regulated
DEG_GutFvsNoGutM_UP <- DEG_GutFvsNoGutM[which(DEG_GutFvsNoGutM$Cutoff%in% "Up"),]
DEG_GutFvsNoGutM_DOWN <- DEG_GutFvsNoGutM[which(DEG_GutFvsNoGutM$Cutoff%in% "Down"),]
DEG_GutFvsNoGutM_UP_annot <-  annotation[which(annotation$gene_id%in%DEG_GutFvsNoGutM_UP$geneID),]
write.csv2(DEG_GutFvsNoGutM_UP_annot, row.names=F, file = paste("path/GutFvsnoGut_DEG_annot_UP.csv") )

#Individual plot figures
#the example is given for GSC/WntA/FoxA and GataB1 from FigureS3 but the list of gene_ID for the other genes is given hereafter 
#list of other genes
gene_list <- read.csv2("path/Gene_Gut_Paper_ID.csv", sep=",")
#Isolate the four target GSC/WntA/FoxA and GataB1
sub_data <- log2.cpm.filtered.norm.df[which(rownames(log2.cpm.filtered.norm.df) %in% c("TRINITY_DN5228_c0_g2","TRINITY_DN12147_c0_g2","TRINITY_DN1633_c3_g1","TRINITY_DN10785_c0_g1")),]
sub_data <- t(as.data.frame(sub_data))
#create a matrix
sub_dat_df <- as.data.frame(matrix(ncol=7,nrow=ncol(sub_data)*nrow(sub_data)))
#named your columns
colnames(sub_dat_df) <- c("target_gene","GeneID","group","sample","log2_cpm_norm","sex", "stage")
sub_dat_df$sample <- rownames(sub_data)

#the follwoing lines work because the samples are ordered the same of you change the order of the sample then update the following
sub_dat_df$group <- targets[which(targets$sample%in%sub_dat_df$sample),"group"]
sub_dat_df$sex <- targets[which(targets$sample%in%sub_dat_df$sample),"sex"]
sub_dat_df$stage <- targets[which(targets$sample%in%sub_dat_df$sample),"stage"]
sub_dat_df$GeneID <- rep(colnames(sub_data),each = nrow(sub_data[,]))
sub_dat_df[which(sub_dat_df$GeneID %in% "TRINITY_DN5228_c0_g2"),"target_gene"] <- "WntA"
sub_dat_df[which(sub_dat_df$GeneID %in% "TRINITY_DN12147_c0_g2"),"target_gene"] <- "FoxA"
sub_dat_df[which(sub_dat_df$GeneID %in% "TRINITY_DN1633_c3_g1"),"target_gene"] <- "GSC"
sub_dat_df[which(sub_dat_df$GeneID %in% "TRINITY_DN10785_c0_g1"),"target_gene"] <- "GataB1"


sub_dat_df$log2_cpm_norm <- sapply(1:nrow(sub_dat_df),function(n){ 
  g <- sub_data[which(rownames(sub_data)%in%sub_dat_df[n,"sample"]),which(colnames(sub_data)%in%sub_dat_df[n,"GeneID"])]
  return(g) })

#Now that you have your data re-arranged for plotting you just want to make sure your values are numeric
sub_dat_df$log2_cpm_norm <-as.numeric(as.character(sub_dat_df$log2_cpm_norm))

#Gut genes

plot <- ggplot(aes(x =factor(stage, levels = c("D0","D1","D2","D3","D4","D5","juvenile", "adult")) , y = log2_cpm_norm, color=factor(sex)), fill=factor(sex), data = sub_dat_df) + 
  facet_wrap(factor(target_gene, levels = c("FoxA","GSC","GataB1","WntA"))~.,ncol=1)+
  geom_point(aes(x =factor(stage, levels = c("D0","D1","D2","D3","D4","D5","juvenile", "adult")), y = log2_cpm_norm), size=3) +
  scale_color_manual(values=c("#E63035","#0C459A","#FECB43"))+
  xlab("Stages")+
  labs(title="Normalized gene expression", subtitle="O.japonicus")+
  theme_light()+theme(legend.position = "none", axis.text.x=element_text(size=10, angle =+30), axis.text.y=element_text(size=12) ,strip.background = element_rect(fill = "white"),strip.text.x.top = element_text(colour ="grey29", face = "bold"))
ggsave(plot = plot, filename= "C:/Users/bqf790/Desktop/R_PostDoc/Plot/Transcriptomic/New_assembly/gut_marker.pdf", dpi = 600, width = 4, height = 12)
ggsave(plot = plot, filename= "C:/Users/bqf790/Desktop/R_PostDoc/Plot/Transcriptomic/New_assembly/gut_marker.tiff", dpi = 600,width = 4, height = 12)
