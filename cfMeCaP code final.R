
################### cfMeCaP manus figures 2024 

library(GenomicRanges)
library(annotatr)
library(dplyr)
library(tidyverse)
library(readxl)
library(survminer)
library(data.table)
library(chromoMap)
library(stringr)
library(GenomicFeatures)
library(ggplot2)
library(RColorBrewer)
library("BSgenome.Hsapiens.UCSC.hg19")
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
library(Hmisc)
library(survival)
library(forestmodel)
library(gghalves)



######  Methylome landscape ###### 

#ctDNA plot
#C1
ggplot(Cohort1_cfDNA_data, aes(x = ID_fac , y = ctDNA, fill=cohort_new)) +
  geom_bar(stat="identity",show.legend = F) +
  scale_fill_manual(values = c("#08519C","#A50F15")) + 
  theme_classic() +
  labs(title= "" ,y="IchorCNA-based ctDNA%", x = "", fill = "") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
        #axis.text.x = element_text(size = 8, colour = "black", angle = 90),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 12, colour = "black"), 
        legend.title = element_text(face="bold"), legend.position = "bottom",
        axis.title.y = element_text(face = "bold", size = 12)) + 
  scale_y_continuous(expand = c(0,0.02), limits = c(0,80)) + 
  facet_grid(~cohort_new, scales = "free_x", space = "free")

#C2
ggplot(Cohort2_comb_cfDNA_data, aes(x = ID_fac , y = ctDNA, fill=cohort_new)) +
  geom_bar(stat="identity",show.legend = F, alpha = 0.9) +
  scale_fill_manual(values = c("#08519C","#FED976", "#FF7F00","#A50F15")) + 
  theme_classic() +
  labs(title= "" ,y="IchorCNA-based ctDNA%", x = "", fill = "") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
        #axis.text.x = element_text(size = 8, colour = "black", angle = 90),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 12, colour = "black"), 
        legend.title = element_text(face="bold"), legend.position = "bottom",
        axis.title.y = element_text(face = "bold", size = 12)) + 
  scale_y_continuous(expand = c(0,0.02), limits = c(0,80)) + 
  facet_grid(~cohort_new, scales = "free_x", space = "free")

###cfDNA plot 
#C1
ggplot(Cohort1_cfDNA_data, aes(x = ID_fac , y = cfDNA, fill=cohort_new)) +
  geom_bar(stat="identity",show.legend = F) +
  scale_fill_manual(values = c("#08519C","#A50F15")) + 
  theme_classic() +
  labs(title= "" ,y="cfDNA concentration (ng/uL)", x = "", fill = "") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
        #axis.text.x = element_text(size = 8, colour = "black", angle = 90),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 12, colour = "black"), 
        legend.title = element_text(face="bold"), legend.position = "bottom",
        axis.title.y = element_text(face = "bold", size = 12)) + 
  scale_y_continuous(expand = c(0,0.02)) + 
  facet_grid(~cohort_new, scales = "free_x", space = "free")

#C2
ggplot(Cohort2_comb_cfDNA_data, aes(x = ID_fac , y = cfDNA, fill=cohort_new)) +
  geom_bar(stat="identity",show.legend = F, alpha = 0.9) +
  scale_fill_manual(values = c("#08519C","#FED976", "#FF7F00","#A50F15")) + 
  theme_classic() +
  labs(title= "" ,y="cfDNA concentration (ng/uL)", x = "", fill = "") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
        #axis.text.x = element_text(size = 8, colour = "black", angle = 90),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 12, colour = "black"), 
        legend.title = element_text(face="bold"), legend.position = "bottom",
        axis.title.y = element_text(face = "bold", size = 12)) + 
  scale_y_continuous(expand = c(0,0.02)) + 
  facet_grid(~cohort_new, scales = "free_x", space = "free")

#####################  Volcano plots #####################  

# define cutoffs for DMR significance 
label_pval_cutoff <- 0.05
label_hyper_cutoff <- 1.5
label_hypo_cutoff <- -1.5

C1A_meth <- readRDS(file = "2022-09-21_TrainControl_highCTDNAmCRPCs_DMRanalysis_allregs_nrpm.rds")

#### VOLCANO PLOT  ####
C1A_meth$region <- paste(C1A_meth$chr, C1A_meth$window_start, sep ="_")

C1A_meth$substring_genes <- str_remove(C1A_meth$genes, pattern = "\\_.*")

C1A_meth$DMR <- "Background"
C1A_meth$DMR[C1A_meth$TvN_log2FC >1.5 & C1A_meth$TvN_adjPval<label_pval_cutoff] <- "HyperDMR"
C1A_meth$DMR[C1A_meth$TvN_log2FC < -1.5 & C1A_meth$TvN_adjPval <label_pval_cutoff] <- "HypoDMR"

C1A_meth$DMR_short <- 0
C1A_meth$DMR_short[C1A_meth$DMR == "HyperDMR"] <- 1
C1A_meth$DMR_short[C1A_meth$DMR == "HypoDMR"] <- -1

C1A_meth$dflabel <- NA
C1A_meth$dflabel[C1A_meth$TvN_log2FC >label_hyper_cutoff & C1A_meth$TvN_adjPval <label_pval_cutoff | C1A_meth$TvN_log2FC < label_hypo_cutoff & C1A_meth$TvN_adjPval <label_pval_cutoff] <- 
  C1A_meth$substring_genes[C1A_meth$TvN_log2FC >label_hyper_cutoff & C1A_meth$TvN_adjPval <label_pval_cutoff | C1A_meth$TvN_log2FC < label_hypo_cutoff & C1A_meth$TvN_adjPval <label_pval_cutoff]

sum(C1A_meth$DMR == "HyperDMR")
# 394
sum(C1A_meth$DMR == "HypoDMR")
# 702
sum(C1A_meth$DMR == "HyperDMR" | C1A_meth$DMR == "HypoDMR")
# 1096

# Volcano plot
C1A_meth %>% 
  ggplot(.,aes(x=TvN_log2FC, y=-log10(TvN_adjPval), col=DMR, label=dflabel)) +
  geom_point(size = 2,  show.legend = F) + #alpha = .2,
  geom_hline(yintercept=-log10(0.05), col="grey", linetype="dashed")+
  geom_vline(xintercept=c(-1.5, 1.5), col="grey", linetype="dashed")+
  xlab("log2FC")+  theme_classic()+  scale_y_continuous(expand = c(0,0), limits = c(0,50)) + 
  ylab("-log10(adjpvalue)")+ theme(axis.title.x = element_text(face = "bold", size = 14), axis.text.x = element_text(face = "bold", size = 14),
                                   axis.title.y = element_text(face = "bold", size = 14),axis.text.y = element_text(face = "bold", size = 14),
                                   panel.background = element_rect(fill='transparent'), #transparent panel bg
                                   plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
                                   panel.grid.major = element_blank(), #remove major gridlines
                                   panel.grid.minor = element_blank(), #remove minor gridlines
                                   legend.background = element_rect(fill='transparent'), #transparent legend bg
                                   legend.box.background = element_rect(fill='transparent'), legend.position = "bottom") + 
  scale_color_manual(values = c("grey", "#CB181D","#2171B5")) 



####################################################################################################################

# Define signficant DMRs for methylome landscape analysis. 
C1A_meth_SigDMRs <- C1A_meth %>%  filter(abs(TvN_log2FC) >1.5 & TvN_adjPval< 0.05) #C1A_meth

C1A_meth_SigDMRs$DMR[C1A_meth_SigDMRs$TvN_log2FC > 1.5 & C1A_meth_SigDMRs$TvN_adjPval< 0.05] <- "HyperDMR"
C1A_meth_SigDMRs$DMR[C1A_meth_SigDMRs$TvN_log2FC < -1.5 & C1A_meth_SigDMRs$TvN_adjPval <0.05] <- "HypoDMR"
C1A_meth_SigDMRs$DMR <- factor(C1A_meth_SigDMRs$DMR, levels = c("HypoDMR", "HyperDMR"))
C1A_meth_SigDMRs$abs_diff <-  abs(C1A_meth_SigDMRs$mCRPC_baseline_nrpm_means - C1A_meth_SigDMRs$control_nrpm_means)


MethSigRegs <- C1A_meth_SigDMRs
MethSigRegs$DMR_short <- ifelse(MethSigRegs$DMR == "Hypermethylated - Cancer enriched", yes = 1, no = 0)
MethSigRegs$text_labs <- NA
MethSigRegs$text_labs  <- MethSigRegs$region

###### LANDSCAPE - DISTRIBUTION OF DMRS: IDEOGRAM ###### 

names(MethSigRegs)
C1A_meth_SigDMRs_ideo_2 <- MethSigRegs[c(1:3, 53)] # Plot dichotomized - hypo vs hyper
names(C1A_meth_SigDMRs_ideo_2) <- c("Chr", "Start", "End", "Value")
C1A_meth_SigDMRs_ideo_2$Chr <- str_remove(string =C1A_meth_SigDMRs_ideo_2$Chr, pattern ="chr")
C1A_meth_SigDMRs_ideo_2$Value <- as.numeric(C1A_meth_SigDMRs_ideo_2$Value)

C1A_meth_SigDMRs_ideo_label <- MethSigRegs[c(1:3, 51)]
names(C1A_meth_SigDMRs_ideo_label) <- c("Chr", "Start", "End", "Value")
C1A_meth_SigDMRs_ideo_label$Chr <- str_remove(string =C1A_meth_SigDMRs_ideo_label$Chr, pattern ="chr")

# CREATE PLOT 
library(RIdeogram)
library(rsvg)

data(human_karyotype, package="RIdeogram")

# SAVE PLOT
ideogram(karyotype = human_karyotype, overlaid = C1A_meth_SigDMRs_ideo_2, label = C1A_meth_SigDMRs_ideo_label, label_type = "polygon", colorset1 = c("#5c95ff","white", "#d62839"))
convertSVG("chromosome.svg", device = "png") 


###### LANDSCAPE - ANNOTATION PLOT: CpG ANNO ###### 

MethSigRegs$CpG_relation_new <- str_remove_all(string = MethSigRegs$CpG_relation, pattern = "hg19_cpg_")

# Clean QSEA annotations 
MethSigRegs$CpG_relation_new[MethSigRegs$CpG_relation_new == "shores, islands" | MethSigRegs$CpG_relation_new == "islands, shores"] <- "islands/shores"
MethSigRegs$CpG_relation_new[MethSigRegs$CpG_relation_new == "shores, shelves" | MethSigRegs$CpG_relation_new == "shelves, shores"] <- "shores/shelves"
MethSigRegs$CpG_relation_new[MethSigRegs$CpG_relation_new == "inter, shelves" | MethSigRegs$CpG_relation_new == "shelves, inter"] <- "shelves/opensea"
MethSigRegs$CpG_relation_new[MethSigRegs$CpG_relation_new == "inter"] <- "opensea"

MethSigRegs$CpG_relation_new_fac <- factor(MethSigRegs$CpG_relation_new, levels = c("islands","islands/shores", "shores","shores/shelves", "shelves","shelves/opensea", "opensea"))

# Set DMR groups 
MethSigRegs$DMR_new <- NA
MethSigRegs$DMR_new[grep(pattern = "Hyper", MethSigRegs$DMR)] <- "HyperDMRs"
MethSigRegs$DMR_new[grep(pattern = "Hypo", MethSigRegs$DMR)] <- "HypoDMRs"

Meth_data <- MethSigRegs %>% dplyr::count(DMR_new, CpG_relation_new_fac) %>%  group_by(DMR_new) %>% 
  mutate(pct=n/sum(n)*100)

ggplot(data = Meth_data, aes(x = DMR_new, fill=CpG_relation_new_fac)) +
  geom_bar(aes(y=pct),show.legend = T, position = "fill", stat = "identity") + 
  labs(y="Percentage of DMRs", x="", fill = "CpG annotation")  + 
  scale_fill_manual(values = c("#08306B", "#08519C", "#2171B5", # "#4292C6", 
                               "#6BAED6" ,"#9ECAE1" ,"#C6DBEF", "#DEEBF7")) + 
  theme_classic() + theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
                          axis.text.x = element_text(size = 12, colour = "black"),
                          axis.text.y = element_text(size = 10, colour = "black"),
                          axis.title.y = element_text(face = "bold", size = 12), legend.position ="bottom") +  
  scale_y_continuous(expand = c(0,0.02), labels=scales::percent)


###### LANDSCAPE - ANNOTATION PLOT: GENE ANNO  ###### 

# Clean QSEA annotations for plotting.
C1A_meth_SigDMRs$genes_relation <- gsub("[^,]*_hg19_", "_hg19_", C1A_meth_SigDMRs$genes)
C1A_meth_SigDMRs$genes_relation <- str_remove_all(string = C1A_meth_SigDMRs$genes_relation, pattern = "_hg19_genes_")
data.frame(C1A_meth_SigDMRs$genes, C1A_meth_SigDMRs$genes_relation)

C1A_meth_SigDMRs$genes_relation_nodups <- gsub("\\b(\\w+)(,\\1\\b)+", "\\1", C1A_meth_SigDMRs$genes_relation, perl = TRUE)
data.frame(C1A_meth_SigDMRs$genes, C1A_meth_SigDMRs$genes_relation, C1A_meth_SigDMRs$genes_relation_nodups)

HyperDMRs <- C1A_meth_SigDMRs %>% filter(TvN_log2FC > 0)
HypoDMRs <- C1A_meth_SigDMRs %>% filter(TvN_log2FC < 0)

No_hyperDMRs <- length(HyperDMRs$chr)
No_hypoDMRs <- length(HypoDMRs$chr)

# Count the number of hyper/hypoDMRs with the specific gene annotation. 
Genes_anno_summary <- rbind(
  data.frame(
    Annotation = c("No annotation" ,"1to5kb", "promoters", "5UTRs", "exons", "introns", "3UTRs" ), 
    Counts = c(length(HyperDMRs$genes_relation[HyperDMRs$genes_relation == ""]),
               length(grep(x = HyperDMRs$genes_relation, pattern = "1to5kb")),
               length(grep(x = HyperDMRs$genes_relation, pattern = "promoters")),
               length(grep(x = HyperDMRs$genes_relation, pattern = "5UTRs")),
               length(grep(x = HyperDMRs$genes_relation, pattern = "exons")),
               length(grep(x = HyperDMRs$genes_relation, pattern = "introns")),
               length(grep(x = HyperDMRs$genes_relation, pattern = "3UTRs"))),
    DMR = c(rep("HyperDMRs", 7))), 
  data.frame(
    Annotation = c("No annotation" ,"1to5kb", "promoters", "5UTRs", "exons", "introns", "3UTRs" ), 
    Counts = c(length(HypoDMRs$genes_relation[HypoDMRs$genes_relation == ""]),
               length(grep(x = HypoDMRs$genes_relation, pattern = "1to5kb")),
               length(grep(x = HypoDMRs$genes_relation, pattern = "promoters")),
               length(grep(x = HypoDMRs$genes_relation, pattern = "5UTRs")),
               length(grep(x = HypoDMRs$genes_relation, pattern = "exons")),
               length(grep(x = HypoDMRs$genes_relation, pattern = "introns")),
               length(grep(x = HypoDMRs$genes_relation, pattern = "3UTRs"))),
    DMR = c(rep("HypoDMRs", 7))))


Genes_anno_summary$perc[Genes_anno_summary$DMR == "HyperDMRs"] <- 
  paste0(round(Genes_anno_summary$Counts[Genes_anno_summary$DMR == "HyperDMRs"]/No_hyperDMRs*100, digits = 1),"%")

Genes_anno_summary$perc[Genes_anno_summary$DMR == "HypoDMRs"] <- 
  paste0(round(Genes_anno_summary$Counts[Genes_anno_summary$DMR == "HypoDMRs"]/No_hypoDMRs*100, digits = 1),"%")

Genes_anno_summary$Annotation <- factor(Genes_anno_summary$Annotation, 
                                        levels = c("No annotation" ,"1to5kb", "promoters", "5UTRs", "exons", "introns", "3UTRs" ))

# color by DMR group
ggplot(Genes_anno_summary, aes(x=Annotation, y=Counts, fill = DMR))+
  geom_bar(stat="identity", color="black", alpha = .8)+ 
  facet_wrap(~DMR, ) + 
  geom_text(aes(label=perc), vjust=-0.25, color="black", size=3)+
  scale_fill_manual(values=c("#A50F15","#08519C"))+  # c("#08519C","#FED976", "#FF7F00","#A50F15")
  theme_classic() + 
  labs(x="", y="No of DMRs") + theme(axis.title.y = element_text(face = "bold"), axis.text.x = element_text(angle = 45)) + 
  scale_y_continuous(expand = c(0,10), limits = c(0,12000))

# color by annotation - this one is used in manus. 
Genes_anno_summary %>% filter(DMR == "HyperDMRs") %>% 
ggplot(., aes(x=Annotation, y=Counts, fill = Annotation))+
  geom_bar(stat="identity", color="black", alpha = 1, show.legend = F)+ 
  #facet_wrap(~DMR, ) + 
  geom_text(aes(label=perc), vjust=-0.25, color="black", size=3)+
  scale_fill_brewer(palette = "Greys") + 
  theme_classic() + 
  labs(x="", y="Number of DMRs ", fill = "Gene annotation", subtitle = "HyperDMRs") + 
  theme(axis.title.y = element_text(face = "bold", size = 12),axis.text.x = element_text(angle = 45, hjust =1), legend.position = "bottom") + 
  scale_y_continuous(expand = c(0,10), limits = c(0,12000))

Genes_anno_summary %>% filter(DMR == "HypoDMRs") %>% 
  ggplot(., aes(x=Annotation, y=Counts, fill = Annotation))+
  geom_bar(stat="identity", color="black", alpha = 1, show.legend = F)+ 
  #facet_wrap(~DMR, ) + 
  geom_text(aes(label=perc), vjust=-0.25, color="black", size=3)+
  scale_fill_brewer(palette = "Greys") + 
  theme_classic() + 
  labs(x="", y="Number of DMRs ", fill = "Gene annotation", subtitle = "HypoDMRs") + 
  theme(axis.title.y = element_text(face = "bold", size = 12),axis.text.x = element_text(angle = 45, hjust =1), legend.position = "bottom") + 
  scale_y_continuous(expand = c(0,10), limits = c(0,12000))

Genes_anno_summary %>% #filter(DMR == "HypoDMRs") %>% 
  ggplot(., aes(x=Annotation, y=Counts, fill = Annotation))+
  geom_bar(stat="identity", color="black", alpha = 1, show.legend = F)+ 
  facet_wrap(~DMR, ) + 
  geom_text(aes(label=perc), vjust=-0.25, color="black", size=3)+
  scale_fill_brewer(palette = "Greys") + 
  theme_classic() + 
  labs(x="", y="Number of DMRs ", fill = "Gene annotation", subtitle = "HypoDMRs") + 
  theme(axis.title.y = element_text(face = "bold", size = 12),axis.text.x = element_text(angle = 45, hjust =1), legend.position = "bottom") + 
  scale_y_continuous(expand = c(0,10), limits = c(0,12000))


####################################################################################################################


######### cfMeCaP signature analyses - by cohort  ################ 

### Load signature regions 
MethSigRegs <- readRDS("2022-12-21_48DMRsSignature_450Ksdfiltered_regionsInfo.rds")

MethSigRegs$region <- rownames(MethSigRegs)

MethSigRegs$region_full <- paste0(MethSigRegs$chr,":", MethSigRegs$window_start, "-", MethSigRegs$window_end)

MethSigRegs$substring_genes <- str_remove(MethSigRegs$genes, pattern = "\\_.*")
MethSigRegs$substring_genes[MethSigRegs$substring_genes == "NA"] <- ""

MethSigRegs$newgene <- MethSigRegs$substring_genes

###############################################################################

##### COHORT 1 - AUH(/HEV)  ###### 

### Read methylation data from patients 
UpdatedCohorts1and2_methData <- readRDS("2022-12-21_48DMRsSignature_450Ksdfiltered_MethLevels.rds")


### Add ctDNA detection
UpdatedCohorts1and2_methData$IchorCNA_detection <- ifelse(test = UpdatedCohorts1and2_methData$ctDNA_fration >= 3, 
                                                          yes = "IchorCNA positive", no = "IchorCNA negative") ## !CHANGE TO 3%

UpdatedCohorts1and2_methData$Group <- factor(UpdatedCohorts1and2_methData$Group, levels = c("NPCC","LPC","HSPC","mCRPC"))

#######  MEAN cfMeCaP METHYLATION  ######

#C1
UpdatedCohorts1and2_methData %>% filter(split == "TRAIN") %>% 
  ggplot(., aes(x = stage , y = DMRs_mean_nrpm)) +
  geom_boxplot(aes(fill = stage),outlier.shape = NA, show.legend = F, alpha = 1) +
  geom_jitter(position=position_jitter(0.3),  size = 3, alpha=1) +
  scale_fill_manual(values = c("#08519C","#A50F15")) +
  theme_classic() + coord_cartesian(ylim = c(0,4))+
  labs(title= "" ,y="mean cfMeCaP methylation (nrpm)", x = "") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
        axis.text.x = element_text(size = 12, colour = "black"),
        axis.text.y = element_text(size = 10, colour = "black"),
        axis.title.y = element_text(face = "bold", size = 12, hjust = .5)) + 
  stat_compare_means(#comparisons = my_comparisons_DMRmean_cohort1, 
    method = "wilcox.test", size = 3) #, label = "p.signif"


library(ggplot2)
library(gghalves)
library(ggpubr)

UpdatedCohorts1and2_methData %>% 
  filter(split == "TRAIN") %>% 
  ggplot(aes(x = stage, y = DMRs_mean_nrpm, fill = stage)) +
  # Half violin on the left
  #geom_half_violin(side = "r", alpha = 1, show.legend = F, scale = "width") +
  geom_violin(alpha = 1, show.legend = F, scale = "width") +
  # Boxplot in the middle
  geom_boxplot(width = 0.15, outlier.shape = NA, show.legend = F, 
               fill = "white", alpha = 1) +
  # Points on the right with jitter
  #geom_half_point(side = "l", color = "black", shape = 19,
  #                size = 1.5, alpha = 1, show.legend = F, 
  #                transformation =position_jitter(width = 0.07, height = 0)) + 
  scale_fill_manual(values = c("#08519C", "#A50F15")) +
  scale_color_manual(values = c("#08519C", "#A50F15")) +
  coord_cartesian(ylim = c(0, 4), clip = "off") +
  labs(y = "mean cfMeCaP methylation (nrpm)", x = "") +
  theme_classic() +
  theme(axis.text.x = element_text(size = 12, colour = "black"),
        axis.text.y = element_text(size = 10, colour = "black"),
        axis.title.y = element_text(face = "bold", size = 12, hjust = 0.5)) +
  stat_compare_means(method = "wilcox.test", size = 3)


my_comparisons_DMRmean_cohort1 <- list(c( "Controls", "RP"),
                                       c("Controls","HSPC"),
                                       c("Controls","mCRPC"))
# C2 

UpdatedCohorts1and2_methData %>% 
  filter(split == "TEST") %>% 
  ggplot(aes(x = stage, y = DMRs_mean_nrpm, fill = stage)) +
  # Half violin on the left
  #geom_half_violin(side = "r", alpha = 1, show.legend = F, scale = "width") +
  geom_violin(alpha = 1, show.legend = F, scale = "width") +
  # Boxplot in the middle
  geom_boxplot(width = 0.15, outlier.shape = NA, show.legend = F, 
               fill = "white", alpha = 1) +
  # Points on the right with jitter
  #geom_half_point(side = "l", color = "black", shape = 19,
  #                size = 1.5, alpha = 1, show.legend = F, 
  #                transformation =position_jitter(width = 0.07, height = 0)) + 
  scale_fill_manual(values = c("#08519C","#FED976", "#FF7F00", "#A50F15", "#67000D")) +
  scale_color_manual(values = c("#08519C","#FED976", "#FF7F00", "#A50F15", "#67000D")) +
  coord_cartesian(ylim = c(0, 4), clip = "off") +
  labs(y = "mean cfMeCaP methylation (nrpm)", x = "") +
  theme_classic() +
  theme(axis.text.x = element_text(size = 12, colour = "black"),
        axis.text.y = element_text(size = 10, colour = "black"),
        axis.title.y = element_text(face = "bold", size = 12, hjust = 0.5)) +
  stat_compare_means(method = "wilcox.test", size = 3)


UpdatedCohorts1and2_methData %>% group_by(split, stage) %>% summarise(cfMeCaP_median=median(DMRs_mean_nrpm)) 


###### Mean methylation of each cfMeCaP region for cohort 1 controls and mCRPC patients #####
mCRPC_col = "#A50F15"
controls_col = "#08519C" 

ggplot(MethSigRegs, aes(x=chr)) + theme_bw() +     
  geom_jitter(aes(y=control_nrpm_means), shape=16, position=position_jitter(0.35),  size = 3, alpha = 1, show.legend = F, col = controls_col) + 
  geom_jitter(aes(y=mCRPC_baseline_nrpm_means), shape=16, position=position_jitter(0.35),  size = 3, alpha = 1, show.legend = F, col = mCRPC_col) +
  theme(axis.text.x = element_blank(), 
        axis.ticks = element_blank(), 
        axis.title.y = element_text(face = "bold", size = 12, vjust = 2),strip.text.x = element_text(size = 12, face = "bold", angle = 0)) + 
  
  scale_y_continuous(expand = c(0,0.2)) + 
  labs(x="", y = "Mean region methylation (nrpm)") + facet_grid(~chr, scales = "free_x", space = "free_x")


### cfMeCaP vs PSA
UpdatedCohorts1and2_methData %>%  filter(split == "TEST" ) %>% #& stage == "mCRPC"
  ggplot(., aes(x=as.factor(mectDNA_detection), # me-ctDNA detection
                y=as.numeric(PSA))) + 
  geom_jitter(aes(), position=position_jitter(0.3), size=3.5, show.legend = TRUE, alpha = 1 ) + #col = factor(M1_mectDNA_detection)
  #scale_color_manual(values = c("#1F78B4", "#E31A1C")) +
  scale_color_manual(values =c("#FED976", "#A50F15")) +
  theme_classic() + scale_fill_brewer(palette = "Greys") + 
  labs(x="",y="serum PSA level", col = "me-ctDNA detection") + 
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
        axis.text.y = element_text(size = 12, colour = "black"),
        axis.text.x = element_text(size = 12, colour = "black"),
        axis.title.y = element_text(face = "bold", size = 12)) + stat_compare_means(method = "wilcox.test", size = 3) + 
  stat_summary(fun = median, geom = "crossbar", width = .5, color = "#525252") + 
  scale_x_discrete(labels=c("0" = "me-ctDNA \n negative", "1" = "me-ctDNA \n positive")) + facet_grid(~Group)

UpdatedCohorts1and2_methData %>%  filter(split == "TRAIN" ) %>% #& stage == "mCRPC"
  ggplot(., aes( x=DMRs_mean_nrpm , y=as.numeric(PSA), col=Group)) + #  methylation
  geom_point(size = 3, show.legend = F, alpha = .8) + # col = factor(M1_mectDNA_detection)
  geom_smooth(method=lm, fullrange=TRUE, color='#2C3E50' ) +
  labs(x="cfMeCaP methylation (nrpm)", y="serum PSA level (ng/mL)") + 
  theme_classic() +  scale_color_manual(values = c("#08519C", "#A50F15")) +
  theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 10, vjust = 0.2), 
        axis.title.x = element_text(face = "bold",size = 10, vjust = 0.5),
        axis.title.y = element_text(face = "bold", size = 10, vjust = 3),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        legend.position = "right", legend.title=element_blank()) +
  stat_cor(method="spearman", cor.coef.name = "rho", label.x.npc = "middle")  

UpdatedCohorts1and2_methData %>%  filter(split == "TEST" ) %>% #& stage == "RP"
  ggplot(., aes( x=DMRs_mean_nrpm , y=as.numeric(PSA), col=Group)) + #  methylation
  geom_point(size = 3, show.legend = F, alpha = .8) + # col = factor(M1_mectDNA_detection)
  geom_smooth(method=lm, fullrange=TRUE, color='#2C3E50' ) +
  labs(x="cfMeCaP methylation (nrpm)", y="serum PSA level (ng/mL)") + 
  theme_classic() + scale_color_manual(values = c("#08519C","#FED976", "#FF7F00", "#A50F15")) +
  theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 10, vjust = 0.2), 
        axis.title.x = element_text(face = "bold",size = 10, vjust = 0.5),
        axis.title.y = element_text(face = "bold", size = 10, vjust = 3),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        legend.position = "right", legend.title=element_blank()) +
  stat_cor(method="spearman", cor.coef.name = "rho", label.x.npc = "middle")  

UpdatedCohorts1and2_methData$metvol_new[UpdatedCohorts1and2_methData$metvol == "0"] <- "M0"
UpdatedCohorts1and2_methData$metvol_new[UpdatedCohorts1and2_methData$metvol == "1"] <- "Bone"
UpdatedCohorts1and2_methData$metvol_new[UpdatedCohorts1and2_methData$metvol == "2"] <- "LN"
UpdatedCohorts1and2_methData$metvol_new[UpdatedCohorts1and2_methData$metvol == "3"] <- "Bone&LN"
UpdatedCohorts1and2_methData$metvol_new[UpdatedCohorts1and2_methData$metvol == "4"] <- "Visceral"

UpdatedCohorts1and2_methData$metvol_new <- factor(UpdatedCohorts1and2_methData$metvol_new, levels=c("M0","Bone","LN","Bone&LN", "Visceral"))

Metvol_comp <- list(c("Bone", "LN"),c("Bone", "Bone&LN"),c("Bone", "Visceral"), 
                    c("LN", "Bone&LN"),c("Visceral", "Bone&LN"), c("LN", "Visceral"))

# Metvol
UpdatedCohorts1and2_methData %>% filter(split == "TRAIN" & Group == "mCRPC") %>% 
  ggplot(., aes(x=metvol_new, # me-ctDNA detection
                y=DMRs_mean_nrpm))+
  geom_boxplot(aes(fill = metvol_new),outlier.shape = NA, show.legend = F, alpha = 1) +
  geom_jitter(position=position_jitter(0.3),  size = 3, alpha=1) +
  theme_classic() +  scale_fill_brewer(palette = "Greys") + 
  labs(title= "" ,y="mean cfMeCaP methylation (nrpm)", x = "") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
        axis.text.x = element_text(size = 9, colour = "black"),
        axis.text.y = element_text(size = 10, colour = "black"),
        axis.title.y = element_text(face = "bold", size = 12)) + 
  #stat_compare_means(comparisons=Metvol_comp,method = "wilcox.test")#label = "p.signif",
  stat_compare_means(method = "kruskal.test", label.x.npc = "middle")

UpdatedCohorts1and2_methData %>% filter(split == "TEST" & Group == "mCRPC") %>% 
  ggplot(., aes(x=metvol_new, # me-ctDNA detection
                y=DMRs_mean_nrpm))+
  geom_boxplot(aes(fill = metvol_new),outlier.shape = NA, show.legend = F, alpha = 1) +
  geom_jitter(position=position_jitter(0.3),  size = 3, alpha=1) +
  theme_classic() +  scale_fill_brewer(palette = "Greys") + 
  labs(title= "" ,y="mean cfMeCaP methylation (nrpm)", x = "") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
        axis.text.x = element_text(size = 9, colour = "black"),
        axis.text.y = element_text(size = 10, colour = "black"),
        axis.title.y = element_text(face = "bold", size = 12)) + 
  stat_compare_means(comparisons=Metvol_comp,method = "wilcox.test")#label = "p.signif",
  #stat_compare_means(method = "kruskal.test", label.x.npc = "middle")

UpdatedCohorts1and2_methData %>% filter(split == "TEST" & Group == "HSPC") %>% 
  ggplot(., aes(x=metvol_new, # me-ctDNA detection
                y=DMRs_mean_nrpm))+
  geom_boxplot(aes(fill = metvol_new),outlier.shape = NA, show.legend = F, alpha = 1) +
  geom_jitter(position=position_jitter(0.3),  size = 3, alpha=1) +
  theme_classic() +  scale_fill_brewer(palette = "Greys") + 
  labs(title= "" ,y="mean cfMeCaP methylation (nrpm)", x = "") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
        axis.text.x = element_text(size = 9, colour = "black"),
        axis.text.y = element_text(size = 10, colour = "black"),
        axis.title.y = element_text(face = "bold", size = 12)) + 
  stat_compare_means(comparisons=list(c("M0", "Bone"), c("M0", "LN"),c("M0", "Bone&LN")),method = "wilcox.test")#label = "p.signif",
#stat_compare_means(method = "kruskal.test", label.x.npc = "middle")


# Barplot . HSPC patients
UpdatedCohorts1and2_methData_HSPC <- UpdatedCohorts1and2_methData %>% filter(Group == "HSPC")

UpdatedCohorts1and2_methData_HSPC$Metvol_grps <- "M1"
UpdatedCohorts1and2_methData_HSPC$Metvol_grps[UpdatedCohorts1and2_methData_HSPC$metvol_new == "M0"] <- "M0"

Barplot_data_gr <- as.data.frame(table(UpdatedCohorts1and2_methData_HSPC$mectDNA_detection, 
                                       UpdatedCohorts1and2_methData_HSPC$Metvol_grps))
totals_gr<- Barplot_data_gr %>%
  group_by(Var2) %>%
  summarise(total=sum(Freq))

#Percentage
ggplot(Barplot_data_gr %>% group_by(Var2) %>%mutate(perc = round(Freq/sum(Freq),2)), aes(x = Var2, y = perc, fill = Var1, cumulative = TRUE)) +
  geom_col(show.legend = FALSE) + theme_bw()+ 
  #scale_fill_manual(values=c("#E31A1C")) + 
  scale_fill_manual(values=c("#1F78B4", "#E31A1C")) + 
  scale_x_discrete(labels=c(paste0(totals_gr$Var2,"\nn=", totals_gr$total)))+
  geom_text(aes(label = paste0(perc*100,"%")), position = position_stack(vjust = 0.5))+
  labs(x = "Phenotype", y="Fraction of patients", fill = "me-ctDNA detection status") + 
  theme(plot.title = element_text(hjust = 0.5, face ="bold"), legend.position = "bottom")


### CONTROLS: Calculate  cutoff again 
# Define function: mean+2*SD
mean_2SD <- function(x){
  return(mean(x)+2*sd(x))
}

# Calculate cutoff for controls (train, n=10)
controls_mean2SD <- subset(Cohort1, Group == "control_train")
Detection_cutoff_controls <- mean_2SD(controls_mean2SD$DMRs_mean_nrpm)


####### cfMeCaP me-ctDNA detection #####
## C1 _______
Cohort1A  <- UpdatedCohorts1and2_methData %>% filter(split == "TRAIN")
Cohort1A$Group2.fac <- factor(Cohort1A$Group, levels = c("NCC", "mCRPC")) 

# BARPLOT 
Barplot_data_gr <- as.data.frame(table(Cohort1A$mectDNA_detection, Cohort1A$Group2.fac))

# ctDNA 
pvalue_gr <- fisher.test(table(Cohort1A$mectDNA_detection,  Cohort1A$Group2.fac), workspace = 2e8)$p
totals_gr<- Barplot_data_gr %>%
  group_by(Var2) %>%
  summarise(total=sum(Freq))

#Percentage
ggplot(Barplot_data_gr %>% group_by(Var2) %>%mutate(perc = round(Freq/sum(Freq),2)), aes(x = Var2, y = perc, fill = Var1, cumulative = TRUE)) +
  geom_col(show.legend = FALSE) + theme_bw()+ scale_fill_manual(values=c("#1F78B4", "#E31A1C")) + 
  scale_x_discrete(labels=c(paste0(totals_gr$Var2,"\nn=", totals_gr$total)))+
  geom_text(aes(label = paste0(perc*100,"%")), position = position_stack(vjust = 0.5))+
  labs(x = "", y="Fraction of patients", fill = "me-ctDNA detection status") + 
  theme(plot.title = element_text(hjust = 0.5, face ="bold"), legend.position = "right")


##  C2 ________
Cohort1B  <- UpdatedCohorts1and2_methData %>% filter(split == "TEST")
Cohort1B$Group2.fac <- factor(Cohort1B$Group, levels = c("NCC","LPC","HSPC", "mCRPC")) 


# BARPLOT 
Barplot_data_gr <- as.data.frame(table(Cohort1B$mectDNA_detection, Cohort1B$Group2.fac))

# ctDNA 
pvalue_gr <- fisher.test(table(Cohort1B$mectDNA_detection,  Cohort1B$Group2.fac), workspace = 2e8)$p
totals_gr<- Barplot_data_gr %>%
  group_by(Var2) %>%
  summarise(total=sum(Freq))

#Percentage
ggplot(Barplot_data_gr %>% group_by(Var2) %>%mutate(perc = round(Freq/sum(Freq),2)), aes(x = Var2, y = perc, fill = Var1, cumulative = TRUE)) +
  geom_col(show.legend = FALSE) + theme_bw()+ scale_fill_manual(values=c("#1F78B4", "#E31A1C")) + 
  scale_x_discrete(labels=c(paste0(totals_gr$Var2,"\nn=", totals_gr$total)))+
  geom_text(aes(label = paste0(perc*100,"%")), position = position_stack(vjust = 0.5))+
  labs(x = "Phenotype", y="Fraction of patients", fill = "me-ctDNA detection status") + 
  theme(plot.title = element_text(hjust = 0.5, face ="bold"), legend.position = "bottom")
# save eps 225* 350


### ICHORCNA DETECTION ####
## C1 _______
Cohort1A  <- UpdatedCohorts1and2_methData %>% filter(split == "TRAIN")
Cohort1A$Group2.fac <- factor(Cohort1A$Group, levels = c("NCC", "mCRPC"))  

# BARPLOT 
Barplot_data_gr <- as.data.frame(table(Cohort1A$IchorCNA_detection, Cohort1A$Group2.fac))

# ctDNA 
pvalue_gr <- fisher.test(table(Cohort1A$IchorCNA_detection,  Cohort1A$Group2.fac), workspace = 2e8)$p
totals_gr<- Barplot_data_gr %>%
  group_by(Var2) %>%
  summarise(total=sum(Freq))

#Percentage
ggplot(Barplot_data_gr %>% group_by(Var2) %>%mutate(perc = round(Freq/sum(Freq),2)), aes(x = Var2, y = perc, fill = Var1, cumulative = TRUE)) +
  geom_col(show.legend = FALSE) + theme_bw()+ scale_fill_manual(values=c("#cfdbd5", "#475b63")) + 
  scale_x_discrete(labels=c(paste0(totals_gr$Var2,"\nn=", totals_gr$total)))+
  geom_text(aes(label = paste0(perc*100,"%")), position = position_stack(vjust = 0.5))+
  labs(x = "", y="Fraction of patients", fill = "me-ctDNA detection status") + 
  theme(plot.title = element_text(hjust = 0.5, face ="bold"), legend.position = "right")
# save eps H225* W200


##  C2 ________
Cohort1B  <- UpdatedCohorts1and2_methData %>% filter(split == "TEST")
Cohort1B$Group2.fac <- factor(Cohort1B$Group, levels = c("NCC","LPC","HSPC", "mCRPC")) 


# BARPLOT 
Barplot_data_gr <- as.data.frame(table(Cohort1B$IchorCNA_detection, Cohort1B$Group2.fac))

# ctDNA 
pvalue_gr <- fisher.test(table(Cohort1B$IchorCNA_detection,  Cohort1B$Group2.fac), workspace = 2e8)$p
totals_gr<- Barplot_data_gr %>%
  group_by(Var2) %>%
  summarise(total=sum(Freq))

#Percentage
ggplot(Barplot_data_gr %>% group_by(Var2) %>%mutate(perc = round(Freq/sum(Freq),2)), aes(x = Var2, y = perc, fill = Var1, cumulative = TRUE)) +
  geom_col(show.legend = FALSE) + theme_bw()+ scale_fill_manual(values=c("#cfdbd5", "#475b63")) + 
  scale_x_discrete(labels=c(paste0(totals_gr$Var2,"\nn=", totals_gr$total)))+
  geom_text(aes(label = paste0(perc*100,"%")), position = position_stack(vjust = 0.5))+
  labs(x = "Phenotype", y="Fraction of patients", fill = "me-ctDNA detection status") + 
  theme(plot.title = element_text(hjust = 0.5, face ="bold"), legend.position = "bottom")
# save eps 225* 350


#### cfMeCaP detection vs IchorCNA detection  ####
Cohort1_PC  <- UpdatedCohorts1and2_methData %>% filter(Group != "NCC")
Cohort1_PC$Group2.fac <- factor(Cohort1_PC$Group, levels = c("LPC","HSPC", "mCRPC")) 

# BARPLOT 
Barplot_data_gr <- as.data.frame(table(Cohort1_PC$mectDNA_detection, Cohort1_PC$IchorCNA_detection))

# ctDNA 
pvalue_gr <- fisher.test(table(Cohort1_PC$mectDNA_detection,  Cohort1_PC$IchorCNA_detection), workspace = 2e8)$p
totals_gr<- Barplot_data_gr %>%
  group_by(Var2) %>%
  summarise(total=sum(Freq))

#Percentage
ggplot(Barplot_data_gr %>% group_by(Var2) %>%mutate(perc = round(Freq/sum(Freq),2)), aes(x = Var2, y = perc, fill = Var1, cumulative = TRUE)) +
  geom_col(show.legend = FALSE) + theme_bw()+ scale_fill_manual(values=c("#1F78B4", "#E31A1C")) + 
  scale_x_discrete(labels=c(paste0(totals_gr$Var2,"\nn=", totals_gr$total)))+
  geom_text(aes(label = paste0(perc*100,"%")), position = position_stack(vjust = 0.5))+
  labs(x = "", y="Fraction of patients", fill = "me-ctDNA detection status") + 
  theme(plot.title = element_text(hjust = 0.5, face ="bold"), legend.position = "bottom")



# Detection overview - cohort 1 mCRPC patients 

detection_cohort <- Cohort1_PC #%>% filter(Cohort == "Cohort1")# & Group  == "LPC"
#table(detection_cohort$mectDNA_detection)
#table(detection_cohort$IchorCNA_detection)

detection_cohort$IchorCNA_detection <- factor(detection_cohort$IchorCNA_detection, levels = c("IchorCNA negative", "IchorCNA positive"))
detection_cohort$mectDNA_detection <- factor(detection_cohort$mectDNA_detection, levels = c(0,1))

detection_table <- table(UpdatedCohorts1and2_methData_HSPC$mectDNA_detection, 
                         UpdatedCohorts1and2_methData_HSPC$Metvol_grps)

detection_table

#library(exact2x2)
mcnemar.exact(detection_table) # if b or c = 0 or if b+c < 25

mcnemar.test(detection_table) 



############  Survival analysis - Cohort 1 mCRPC patients, n=48 ############  
library(forestmodel)
library(survival)
library(survminer)


C1and2_SurvData$cfMeCaP_methylation_GR<- cut(C1and2_SurvData$cfMeCaP_methylation, breaks=c(-Inf,median(C1and2_SurvData$cfMeCaP_methylation),+Inf), labels=c("low","high")) 
C1and2_SurvData$IchorCNA_ctDNAfrac_GR<- cut(C1and2_SurvData$ctDNA_fraction, breaks=c(-Inf,median(C1and2_SurvData$ctDNA_fraction),+Inf), labels=c("low","high")) 
C1and2_SurvData$Baseline_PSA_GR <- cut(C1and2_SurvData$Baseline_PSA , breaks=c(-Inf,median(C1and2_SurvData$Baseline_PSA ),+Inf), labels=c(0,1)) 


##### 
# ___Cox - univariate - PFS ___
res.cox_pfs <-coxph(Surv(PFS_mo,PFS_status)~ctDNA_fraction, data=C1and2_SurvData)

cox.zph(res.cox_pfs) # Test for proportional hazards. 
summary(res.cox_pfs)


# plot martingale residuals
res <- residuals(res.cox_pfs, type = "martingale")

plot(C1and2_SurvData$cfMeCaP_methylation, res, xlab = "X", ylab = "Martingale residuals")
lines(lowess(C1and2_SurvData$cfMeCaP_methylation, res), col = "red", lwd = 2)
abline(h = 0, lty = 2)


# ___ Kaplan Meier - PFS ___
pfs_frag<-survfit(Surv(PFS_mo,PFS_status)~cfMeCaP_detection, data=C1and2_SurvData)

ggsurvplot(pfs_frag, data=C1and2_SurvData,
           combine=TRUE, 
           risk.table=TRUE, 
           conf.int=FALSE, 
           conf.int.style="ribbon",                     # step or ribbon
           censor=TRUE,                                 #TRUE: show censoring (length of FU)
           tables.theme=theme_cleantable(),
           title = "",
           legend.title = "",legend.labs = c("me-ctDNA negative", "me-ctDNA positive"), #legend.labs = c("Low cfMeCaP methylation", "High cfMeCaP methylation"),  c("2-30% ctDNA", "30-100% ctDNA")
           xlab="Progression-free survival (months)",
           ylab="PFS probability",
           palette=c("dodgerblue3","red2", "darkgreen", "purple", "black", "orange", "pink"),                                # other: Set1, jco,
           pval.method = TRUE,
           pval=TRUE, 
           pval.coord = c(50,0.9), #xlim = c(0,40),  
           pval.method.coord = c(50,1), 
           surv.median.line = "hv")  



#OS 
res.cox_os <-coxph(Surv(OS_mo, OS_status)~ctDNA_fraction, data=C1and2_SurvData)

cox.zph(res.cox_os) # Test for proportional hazards. 
summary(res.cox_os)

os_frag<-survfit(Surv(OS_mo, OS_status)~cfMeCaP_methylation_GR, data=C1and2_SurvData)

ggsurvplot(os_frag, data=C1and2_SurvData,
           combine=TRUE, 
           risk.table=TRUE, 
           conf.int=FALSE, 
           conf.int.style="ribbon",                     # step or ribbon
           censor=TRUE,                                 #TRUE: show censoring (length of FU)
           tables.theme=theme_cleantable(),
           title = "",
           legend.title = "",legend.labs = c("Low", "High"), #legend.labs = c("Low cfMeCaP methylation", "High cfMeCaP methylation"),
           xlab="Overall survival (months)",
           ylab="OS probability",
           palette=c("dodgerblue3","red2", "darkgreen", "purple", "black", "orange", "pink"),                                # other: Set1, jco,
           pval.method = TRUE,
           pval=TRUE, 
           pval.coord = c(30,0.9), #xlim = c(0,40),  
           pval.method.coord = c(30,1), 
           surv.median.line = "hv")  



#### ___Forest plot - univariate cox - PFS ___#### 
Cohort1_SurvData_forest <- C1and2_SurvData %>% 
  #filter(treatment == "Enza") %>% 
  transmute(PFS_mo, PFS_status, 
            radioPFS_mo, radioPFS_status, 
            OS_mo, OS_status,
            cfMeCaP = cfMeCaP_methylation, 
            IchorCNA = ctDNA_fraction,
            IchorCNA_2groups = IchorCNA_ctDNAfrac_2groups, IchorCNA_3groups = IchorCNA_ctDNAfrac_3groups, 
            PSA = Baseline_PSA,  
            Mets = factor(Met_vol,labels = c("Bone", "LN", "Bone/LN", "Visceral")), 

            cfMeCaP_dichotomized =factor(cfMeCaP_methylation_GR, labels = c("Low cfMeCaP methylation", "High cfMeCaP methylation")),
            IchorCNA_ctDNA_dichotomized =factor(IchorCNA_ctDNAfrac_GR, labels = c("Low ctDNA%", "High ctDNA%")))



# Combine bone and lymph node groups
Cohort1_SurvData_forest$Metvol_new[Cohort1_SurvData_forest$Mets %in% c("Bone", "Lymph node", "Bone and lymph node")] <-  "Bone and/or lymph node"
Cohort1_SurvData_forest$Metvol_new[Cohort1_SurvData_forest$Mets == "Visceral"] <-  "Visceral"                          
Cohort1_SurvData_forest$Metvol_new[Cohort1_SurvData_forest$Mets == "Unknown"] <-  "Unknown"                          


vars_for_table <-c("cfMeCaP", "IchorCNA", "PSA", "Mets")

# define endpoints 
univ_formulas <- sapply(vars_for_table, function(x) as.formula(paste('Surv(PFS_mo, PFS_status)~', x)))

#making a list of models
univ_models <- lapply(univ_formulas, function(x){coxph(x,data=Cohort1_SurvData_forest)})

#extract data 
univ_results <- lapply(univ_models,function(x){return(exp(cbind(coef(x),confint(x))))})

# print plot
print(forest_model(model_list = univ_models,covariates = vars_for_table,merge_models =T))


# test assumptions for each variable!
res.cox_pfs <-coxph(Surv(PFS_mo, PFS_status)~cfMeCaP + PSA + Mets, data=Cohort1_SurvData_forest)

cox.zph(res.cox_pfs)  
summary(res.cox_pfs)

# adjust p-values 
p.adjust(p = c(), method="BH")


#___Forest plot - Multivariate cox - PFS 

# Select significant (after BH adj.) variables  in univariate analysis
print(forest_model(coxph(Surv(PFS_mo, PFS_status)~cfMeCaP + PSA + Mets, data=Cohort1_SurvData_forest)))


#### ___Forest plot - univariate cox - OS ___#### 

# Test proportional hazards assumptions 
res.cox_os <-coxph(Surv(OS_mo, OS_status)~IchorCNA, data=Cohort1_SurvData_forest)

cox.zph(res.cox_os) 
summary(res.cox_os)


### Create forest plot

vars_for_table_OS <-c("cfMeCaP", "IchorCNA", "PSA", "Mets")
univ_formulas_OS <- sapply(vars_for_table_OS, function(x) as.formula(paste('Surv(OS_mo, OS_status)~', x)))
univ_models_OS <- lapply(univ_formulas_OS, function(x){coxph(x,data=Cohort1_SurvData_forest)})
print(forest_model(model_list = univ_models_OS,covariates = vars_for_table_OS,merge_models =T))

# BH adjust p-vals 
p.adjust(p = c(0.000356, , 0.0617, 0.9437, 0.8834, 0.0033, 0.9612), method="BH")


#### ___Forest plot - Multivariate cox - OS ___#### 
print(forest_model(coxph(Surv(OS_mo, OS_status)~mean_cfMeCaP_methylation + 
                           PSA_level + 
                           Metastatic_volume, data=Cohort1_SurvData_forest)))



####################################################################################################################

############################ COHORT 3 - Chen et al.############################ 

# Load cfMeCaP regions 
#Cohort2_cfMeCaPregions <- readRDS("E:/molpros/faststorage/BACKUP/Karina/N282_cfMeDIP-seq/Scripts/QSEA/2023_EGAS00001005522run/2023-06-13_Cohort3_cfMeCaPregions_qsea300bp.rds")
Cohort3_cfMeCaPregions <- readRDS("2023-06-13_Cohort2_cfMeCaPregions_qsea300bp.rds")

#Load Patient info to identify matched samples 
Cohort3_overview <- read_excel("Cohort3_overview.xlsx")

# PREPARE DATA - analyses will only be based on RP samples and mCRPC baseline samples. 

# Combine data 
Cohort3_cfMeCaPregions_info <- cbind(Cohort3_cfMeCaPregions, 
                                     study = Cohort3_overview$Study[match(rownames(Cohort3_cfMeCaPregions), Cohort3_overview$InHouseID_2)], 
                                     TimePoint = Cohort3_overview$Timepoint[match(rownames(Cohort3_cfMeCaPregions), Cohort3_overview$InHouseID_2)],
                                     MeDIP_sample = Cohort3_overview$MeDIPsampleID[match(rownames(Cohort3_cfMeCaPregions), Cohort3_overview$InHouseID_2)],
                                     Patient_TRUE = Cohort3_overview$Patient[match(rownames(Cohort3_cfMeCaPregions), Cohort3_overview$InHouseID_2)])   


# Filter samples
Cohort3_cfMeCaPregions_info_RP_mCRPCbase <- subset(Cohort3_cfMeCaPregions_info, TimePoint %in% c("lPC", "mCRPC-base", "Baseline"))

Cohort3_cfMeCaPregions_info_RP_mCRPCbase$timepoint_new <- ifelse(Cohort3_cfMeCaPregions_info_RP_mCRPCbase$TimePoint == "lPC", yes = "LPC", no = "mCRPC")
Cohort3_cfMeCaPregions_info_RP_mCRPCbase$timepoint_new <- factor(Cohort3_cfMeCaPregions_info_RP_mCRPCbase$timepoint_new, levels = c("LPC", "mCRPC"))


# Calculate mean cfMeCaP48 methylation per patient
Cohort3_cfMeCaPregions_info_RP_mCRPCbase$DMRs_mean_nrpm <- apply(Cohort3_cfMeCaPregions_info_RP_mCRPCbase[1:48], 1, mean, na.rm=T)


# ctDNA detection: Test for ctDNA detection in the remaining samples 
ctDNA_detection_score <- c()

for (i in 1:length(rownames(Cohort3_cfMeCaPregions_info_RP_mCRPCbase))) {
  DMR_score <- ifelse(test = Cohort3_cfMeCaPregions_info_RP_mCRPCbase[i,"DMRs_mean_nrpm"] > Detection_cutoff_controls  ,yes = 1, no = 0)
  ctDNA_detection_score <- c(ctDNA_detection_score,DMR_score)
}

Cohort3_cfMeCaPregions_info_RP_mCRPCbase$mectDNA_detection <- ctDNA_detection_score

##### Figure 2C: Plot mean cfMeCaP methylation ##### 
ggplot(Cohort3_cfMeCaPregions_info_RP_mCRPCbase, aes(x = timepoint_new , y = DMRs_mean_nrpm)) +
  geom_boxplot(aes(fill = timepoint_new),outlier.shape = NA, show.legend = F, alpha = 1) +
  geom_jitter(position=position_jitter(0.3),  size = 3, alpha=1) +
  scale_fill_manual(values = c("#FED976", "#A50F15")) +
  theme_classic() + 
  labs(title= "" ,y="mean cfMeCaP methylation (nrpm)", x = "") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
        axis.text.x = element_text(size = 9, colour = "black"),
        axis.text.y = element_text(size = 10, colour = "black"),
        axis.title.y = element_text(face = "bold", size = 12)) + 
  stat_compare_means( method = "wilcox.test") #label = "p.signif",

# violin
Cohort3_cfMeCaPregions_info_RP_mCRPCbase %>% 
  ggplot(aes(x = timepoint_new, y = DMRs_mean_nrpm, fill = timepoint_new)) +
  # Half violin on the left
  #geom_half_violin(side = "r", alpha = 1, show.legend = F, scale = "width") +
  geom_violin(alpha = 1, show.legend = F, scale = "width") +
  # Boxplot in the middle
  geom_boxplot(width = 0.15, outlier.shape = NA, show.legend = F, 
               fill = "white", alpha = 1) +
  # Points on the right with jitter
  #geom_half_point(side = "l", color = "black", shape = 19,
  #                size = 1.5, alpha = 1, show.legend = F, 
  #                transformation =position_jitter(width = 0.07, height = 0)) + 
  scale_fill_manual(values = c("#FED976", "#A50F15")) +
  scale_color_manual(values = c("#FED976", "#A50F15")) +
  coord_cartesian(ylim = c(0, 4), clip = "off") +
  labs(y = "mean cfMeCaP methylation (nrpm)", x = "") +
  theme_classic() +
  theme(axis.text.x = element_text(size = 12, colour = "black"),
        axis.text.y = element_text(size = 10, colour = "black"),
        axis.title.y = element_text(face = "bold", size = 12, hjust = 0.5)) +
  stat_compare_means(method = "wilcox.test", size = 3)


# save eps H225 x W240

Cohort3_cfMeCaPregions_info_RP_mCRPCbase %>% 
  group_by(timepoint_new) %>% 
  summarise(median=median(DMRs_mean_nrpm)) 

##### Figure 3C: cfMeCaP me-ctDNA detection ##### 

# BARPLOT 
Barplot_data_gr <- as.data.frame(table(Cohort3_cfMeCaPregions_info_RP_mCRPCbase$mectDNA_detection, Cohort3_cfMeCaPregions_info_RP_mCRPCbase$timepoint_new))

# ctDNA 
pvalue_gr <- fisher.test(table(Cohort3_cfMeCaPregions_info_RP_mCRPCbase$mectDNA_detection, Cohort3_cfMeCaPregions_info_RP_mCRPCbase$timepoint_new), workspace = 2e8)$p
totals_gr<- Barplot_data_gr %>%
  group_by(Var2) %>%
  summarise(total=sum(Freq))

#Percentage
ggplot(Barplot_data_gr %>% group_by(Var2) %>%mutate(perc = round(Freq/sum(Freq),2)), aes(x = Var2, y = perc, fill = Var1, cumulative = TRUE)) +
  geom_col(show.legend = FALSE) + theme_bw()+ scale_fill_manual(values=c("#1F78B4", "#E31A1C")) + 
  scale_x_discrete(labels=c(paste0(totals_gr$Var2,"\nn=", totals_gr$total)))+
  geom_text(aes(label = paste0(perc*100,"%")), position = position_stack(vjust = 0.5))+
  labs(x = "", y="Fraction of patients", fill = "me-ctDNA detection status") + 
  theme(plot.title = element_text(hjust = 0.5, face ="bold"), legend.position = "bottom")


### me-ctDNA vs clinical parameters
C3_mCRPC_info <- as.data.frame(read_excel("Cohort 2_Sampleinfo_clinical_updated.xlsx", sheet = "Full_initial_clinData"))

Cohort2_cfMeCaPregions_info_RP_mCRPCbase_IniClinData <- Cohort2_cfMeCaPregions_info_RP_mCRPCbase %>% 
  left_join(., C3_mCRPC_info, by=c("MeDIP_sample" = "MeDIPsampleID") )


Cohort2_cfMeCaPregions_info_RP_mCRPCbase_IniClinData$Pre_Treatment_PSA <- as.numeric(str_replace(Cohort2_cfMeCaPregions_info_RP_mCRPCbase_IniClinData$Pre_Treatment_PSA, 
                                                                                                 pattern = ",", replacement = "."))

Cohort2_cfMeCaPregions_info_RP_mCRPCbase_IniClinData$Metvol <- "Unknown"
Cohort2_cfMeCaPregions_info_RP_mCRPCbase_IniClinData$Metvol[Cohort2_cfMeCaPregions_info_RP_mCRPCbase_IniClinData$Bone_mets == "Yes"] <- "Bone"
Cohort2_cfMeCaPregions_info_RP_mCRPCbase_IniClinData$Metvol[Cohort2_cfMeCaPregions_info_RP_mCRPCbase_IniClinData$Lung_mets == "Yes"& Cohort2_cfMeCaPregions_info_RP_mCRPCbase_IniClinData$Liver_mets == "No"] <- "Lung"
Cohort2_cfMeCaPregions_info_RP_mCRPCbase_IniClinData$Metvol[Cohort2_cfMeCaPregions_info_RP_mCRPCbase_IniClinData$Lung_mets == "No"& Cohort2_cfMeCaPregions_info_RP_mCRPCbase_IniClinData$Liver_mets == "Yes"] <- "Liver"
Cohort2_cfMeCaPregions_info_RP_mCRPCbase_IniClinData$Metvol[Cohort2_cfMeCaPregions_info_RP_mCRPCbase_IniClinData$Lung_mets == "Yes"& Cohort2_cfMeCaPregions_info_RP_mCRPCbase_IniClinData$Liver_mets == "Yes"] <- "Lung & Liver"


# PSA vs cfMeCaP
Cohort2_cfMeCaPregions_info_RP_mCRPCbase_IniClinData %>% 
  ggplot(., aes(x=as.factor(mectDNA_detection), # me-ctDNA detection
                y=as.numeric(Pre_Treatment_PSA))) + 
  geom_jitter(aes(), position=position_jitter(0.3), size=3.5, show.legend = TRUE, alpha = 1 ) + #col = factor(M1_mectDNA_detection)
  #scale_color_manual(values = c("#1F78B4", "#E31A1C")) +
  scale_color_manual(values =c("#FED976", "#A50F15")) +
  theme_classic() + scale_fill_brewer(palette = "Greys") + 
  labs(x="",y="serum PSA level", col = "me-ctDNA detection") + 
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
        axis.text.y = element_text(size = 12, colour = "black"),
        axis.text.x = element_text(size = 12, colour = "black"),
        axis.title.y = element_text(face = "bold", size = 12)) + stat_compare_means(method = "wilcox.test", size = 3) + 
  stat_summary(fun = median, geom = "crossbar", width = .5, color = "#525252") + 
  scale_x_discrete(labels=c("0" = "me-ctDNA \n negative", "1" = "me-ctDNA \n positive")) + facet_grid(~timepoint_new)

Cohort2_cfMeCaPregions_info_RP_mCRPCbase_IniClinData %>%  filter(timepoint_new == "mCRPC") %>% 
  ggplot(., aes( x=DMRs_mean_nrpm , y=as.numeric(Pre_Treatment_PSA))) + #  methylation
  geom_point(size = 3, show.legend = F, alpha = .8) + # col = factor(M1_mectDNA_detection)
  geom_smooth(method=lm, fullrange=TRUE, color='#2C3E50' ) +
  labs(x="cfMeCaP methylation (nrpm)", y="serum PSA level (ng/mL)") + 
  theme_classic() + 
  theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 10, vjust = 0.2), 
        axis.title.x = element_text(face = "bold",size = 10, vjust = 0.5),
        axis.title.y = element_text(face = "bold", size = 10, vjust = 3),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        legend.position = "right", legend.title=element_blank()) +
  stat_cor(method="spearman", cor.coef.name = "rho", label.x.npc = "middle")  

Metvol_comp <- list(c("Bone", "Liver"),c("Bone", "Lung"),c("Bone", "Lung & Liver"))

# Metvol
Cohort2_cfMeCaPregions_info_RP_mCRPCbase_IniClinData %>% filter(timepoint_new == "mCRPC") %>% 
  ggplot(., aes(x=Metvol, # me-ctDNA detection
                y=DMRs_mean_nrpm))+
  geom_boxplot(aes(fill = Metvol),outlier.shape = NA, show.legend = F, alpha = 1) +
  geom_jitter(position=position_jitter(0.3),  size = 3, alpha=1) +
  theme_classic() +  scale_fill_brewer(palette = "Greys") + 
  labs(title= "" ,y="mean cfMeCaP methylation (nrpm)", x = "") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
        axis.text.x = element_text(size = 9, colour = "black"),
        axis.text.y = element_text(size = 10, colour = "black"),
        axis.title.y = element_text(face = "bold", size = 12)) + 
  stat_compare_means(comparisons=Metvol_comp,method = "wilcox.test") #label = "p.signif",


Barplot_data_gr <- as.data.frame(table(Cohort2_cfMeCaPregions_info_RP_mCRPCbase_IniClinData$mectDNA_detection[Cohort2_cfMeCaPregions_info_RP_mCRPCbase_IniClinData$timepoint_new == "mCRPC"], 
                                       Cohort2_cfMeCaPregions_info_RP_mCRPCbase_IniClinData$Metvol[Cohort2_cfMeCaPregions_info_RP_mCRPCbase_IniClinData$timepoint_new == "mCRPC"]))
totals_gr<- Barplot_data_gr %>%
  group_by(Var2) %>%
  summarise(total=sum(Freq))

#Percentage
ggplot(Barplot_data_gr %>% group_by(Var2) %>%mutate(perc = round(Freq/sum(Freq),2)), aes(x = Var2, y = perc, fill = Var1, cumulative = TRUE)) +
  geom_col(show.legend = FALSE) + theme_bw()+ scale_fill_manual(values=c("#1F78B4", "#E31A1C")) + 
  scale_x_discrete(labels=c(paste0(totals_gr$Var2,"\nn=", totals_gr$total)))+
  geom_text(aes(label = paste0(perc*100,"%")), position = position_stack(vjust = 0.5))+
  labs(x = "Phenotype", y="Fraction of patients", fill = "me-ctDNA detection status") + 
  theme(plot.title = element_text(hjust = 0.5, face ="bold"), legend.position = "bottom")


##### SURVIVAL ANALYSES - CRPC PATIENTS
# Prepare clinical FU data. In the full data set, mCRPC patients from 4 different substudies are in included (VPC, VPC-V, WCDT and Barrier).
# Only patients from the VPC and VPC-V study have the necessary FU available and will be included in the survival analyses. 
# Mutation-based ctDNA% are available for the same patients (n=72)

library(readxl)

### VPC-V study
C3_VPCV_info <- as.data.frame(read_excel("Cohort 2_Sampleinfo_clinical_updated.xlsx", sheet = "VPC-V"))
C3_VPCV_methdata <- Cohort2_cfMeCaPregions_info_RP_mCRPCbase %>% filter(timepoint_new == "mCRPC" & MeDIP_sample %in% C3_VPCV_info$MeDIPsampleID) 

C3_VPCV_methdata_clinFU <- cbind(C3_VPCV_methdata, 
                                 ctDNA =     C3_VPCV_info$ctDNA_frac[match(C3_VPCV_methdata$MeDIP_sample, C3_VPCV_info$MeDIPsampleID)],
                                 PSA =       C3_VPCV_info$PSA[match(C3_VPCV_methdata$MeDIP_sample, C3_VPCV_info$MeDIPsampleID)],
                                 Bone_mets = C3_VPCV_info$`Bone metastases`[match(C3_VPCV_methdata$MeDIP_sample, C3_VPCV_info$MeDIPsampleID)],
                                 Lung_mets = C3_VPCV_info$`Lung metastases`[match(C3_VPCV_methdata$MeDIP_sample, C3_VPCV_info$MeDIPsampleID)],
                                 Liver_mets = C3_VPCV_info$`Liver metastases`[match(C3_VPCV_methdata$MeDIP_sample, C3_VPCV_info$MeDIPsampleID)],
                                 LDH = C3_VPCV_info$`LDH / ULN`[match(C3_VPCV_methdata$MeDIP_sample, C3_VPCV_info$MeDIPsampleID)],
                                 ALP = C3_VPCV_info$`ALP / ULN`[match(C3_VPCV_methdata$MeDIP_sample, C3_VPCV_info$MeDIPsampleID)],
                                 prog_psa =   C3_VPCV_info$Pgr_days[match(C3_VPCV_methdata$MeDIP_sample, C3_VPCV_info$MeDIPsampleID)],
                                 prog_censored = C3_VPCV_info$Pgr_Censored[match(C3_VPCV_methdata$MeDIP_sample, C3_VPCV_info$MeDIPsampleID)],
                                 os_days =         C3_VPCV_info$OS_days[match(C3_VPCV_methdata$MeDIP_sample, C3_VPCV_info$MeDIPsampleID)],
                                 os_censored =     C3_VPCV_info$OS_Censored[match(C3_VPCV_methdata$MeDIP_sample, C3_VPCV_info$MeDIPsampleID)],
                                 Best_PSA_response = C3_VPCV_info$`Best response within first 12 weeks (% change)`[match(C3_VPCV_methdata$MeDIP_sample, C3_VPCV_info$MeDIPsampleID)])


C3_VPCV_methdata_clinFU$prog_psa_mo <- C3_VPCV_methdata_clinFU$prog_psa/30.5 
C3_VPCV_methdata_clinFU$os_mo <- C3_VPCV_methdata_clinFU$os_days/30.5 

# add pgr/os status columns 
C3_VPCV_methdata_clinFU$prog_status <- ifelse(C3_VPCV_methdata_clinFU$prog_censored == "Yes", yes = 0, no = 1)
C3_VPCV_methdata_clinFU$os_status <- ifelse(C3_VPCV_methdata_clinFU$os_censored == "Yes", yes = 0, no = 1)


### VPC
#setwd("U:/Kandidat/MOMA/cfMeDIP-seq paper 2022/Data/")
C3_VPC_info <- as.data.frame(read_excel("Cohort 2_Sampleinfo_clinical_updated.xlsx", sheet = "VPC"))
colnames(C3_VPC_info)[colnames(C3_VPC_info) %in% c("Published_Censored...27", "Published_Censored...29")] <- c("Progression_status", "Death_status")

# Subset VPC-set from methylation data 
C3_VPC_methdata <- Cohort2_cfMeCaPregions_info_RP_mCRPCbase %>% filter(timepoint_new == "mCRPC" & MeDIP_sample %in% C3_VPC_info$MeDIPsampleID) 

C3_VPC_methdata_clinFU <- cbind(C3_VPC_methdata, 
                                ctDNA =     C3_VPC_info$ctDNA_frac[match(C3_VPC_methdata$MeDIP_sample, C3_VPC_info$MeDIPsampleID)],
                                PSA =       C3_VPC_info$PSA[match(C3_VPC_methdata$MeDIP_sample, C3_VPC_info$MeDIPsampleID)],
                                Bone_mets = C3_VPC_info$`Bone metastases`[match(C3_VPC_methdata$MeDIP_sample, C3_VPC_info$MeDIPsampleID)],
                                Lung_mets = C3_VPC_info$`Lung metastases`[match(C3_VPC_methdata$MeDIP_sample, C3_VPC_info$MeDIPsampleID)],
                                Liver_mets = C3_VPC_info$`Liver metastases`[match(C3_VPC_methdata$MeDIP_sample, C3_VPC_info$MeDIPsampleID)],
                                LDH = C3_VPC_info$`LDH / ULN`[match(C3_VPC_methdata$MeDIP_sample, C3_VPC_info$MeDIPsampleID)],
                                ALP = C3_VPC_info$`ALP / ULN`[match(C3_VPC_methdata$MeDIP_sample, C3_VPC_info$MeDIPsampleID)],
                                prog_psa =   C3_VPC_info$Pgr_days[match(C3_VPC_methdata$MeDIP_sample, C3_VPC_info$MeDIPsampleID)],
                                prog_censored = C3_VPC_info$Pgr_Censored[match(C3_VPC_methdata$MeDIP_sample, C3_VPC_info$MeDIPsampleID)],
                                os_days =         C3_VPC_info$OS_days[match(C3_VPC_methdata$MeDIP_sample, C3_VPC_info$MeDIPsampleID)],
                                os_censored =     C3_VPC_info$OS_Censored[match(C3_VPC_methdata$MeDIP_sample, C3_VPC_info$MeDIPsampleID)],
                                Best_PSA_response = C3_VPC_info$`Best response within first 12 weeks (% change)`[match(C3_VPC_methdata$MeDIP_sample, C3_VPC_info$MeDIPsampleID)])


C3_VPC_methdata_clinFU$prog_psa_mo <- C3_VPC_methdata_clinFU$prog_psa/30.5 
C3_VPC_methdata_clinFU$os_mo <- C3_VPC_methdata_clinFU$os_days/30.5 

# add pgr/os status columns 
C3_VPC_methdata_clinFU$prog_status <- ifelse(C3_VPC_methdata_clinFU$prog_censored == "Yes", yes = 0, no = 1)
C3_VPC_methdata_clinFU$os_status <- ifelse(C3_VPC_methdata_clinFU$os_censored == "Yes", yes = 0, no = 1)

### combine VPC and VPC-V data 
C3_mCRPC_methdata_clinFU <- rbind(C3_VPCV_methdata_clinFU, C3_VPC_methdata_clinFU)


##### Figure S6C: Correlate ctDNA% and me-ctDNA levels ##### 

ggplot(C3_mCRPC_methdata_clinFU, aes(x=ctDNA  , y=DMRs_mean_nrpm )) + 
  geom_point(size = 2, show.legend = T) +
  #scale_color_manual(values = c("#FED976", "#FF7F00", "#67000D")) +
  geom_smooth(method=lm, fullrange=TRUE, color='#2C3E50' ) +
  labs(x = "ctDNA%", y = "cfMeCaP48 methylation (nrpm)", title = "") + 
  theme_classic() + 
  theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 10, vjust = 0.2), 
        axis.title.x = element_text(face = "bold",size = 10, vjust = 0.5),
        axis.title.y = element_text(face = "bold", size = 10, vjust = 3),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        legend.position = "right", legend.title=element_blank()) +
  stat_cor(method="spearman", cor.coef.name = "rho", label.x = 1, label.y = 4, size = 4)


##### Figure S6C: Metvol ##### 

C3_mCRPC_methdata_clinFU$Metvol <- "Unknown"
C3_mCRPC_methdata_clinFU$Metvol[C3_mCRPC_methdata_clinFU$Bone_mets == "Yes"] <- "Bone"
C3_mCRPC_methdata_clinFU$Metvol[C3_mCRPC_methdata_clinFU$Lung_mets == "Yes"& C3_mCRPC_methdata_clinFU$Liver_mets == "No"] <- "Lung"
C3_mCRPC_methdata_clinFU$Metvol[C3_mCRPC_methdata_clinFU$Lung_mets == "No"& C3_mCRPC_methdata_clinFU$Liver_mets == "Yes"] <- "Liver"
C3_mCRPC_methdata_clinFU$Metvol[C3_mCRPC_methdata_clinFU$Lung_mets == "Yes"& C3_mCRPC_methdata_clinFU$Liver_mets == "Yes"] <- "Lung & Liver"


##### Figure 3F:  Mutations-based ctDNA detection per group ##### 
C3_mCRPC_methdata_clinFU$ctDNA_detection <- ifelse(test = C3_mCRPC_methdata_clinFU$ctDNA > 0, yes = "ctDNA positive", no = "ctDNA negative")

C3_mCRPC_methdata_clinFU$mectDNA_detection <- ifelse(test = C3_mCRPC_methdata_clinFU$DMRs_mean_nrpm > Detection_cutoff_controls, yes = "mectDNA positive", no = "ctDNA negative")
  

C3_mCRPC_methdata_clinFU$group <- rep("mCRPC", 72)

Barplot_data_gr <- as.data.frame(table(C3_mCRPC_methdata_clinFU$ctDNA_detection, C3_mCRPC_methdata_clinFU$group))
  
# ctDNA 
pvalue_gr <- fisher.test(table(C3_mCRPC_methdata_clinFU$ctDNA_detection, C3_mCRPC_methdata_clinFU$group), workspace = 2e8)$p
totals_gr<- Barplot_data_gr %>%
  group_by(Var2) %>%
  summarise(total=sum(Freq))

#Percentage
ggplot(Barplot_data_gr %>% group_by(Var2) %>%mutate(perc = round(Freq/sum(Freq),2)), aes(x = Var2, y = perc, fill = Var1, cumulative = TRUE)) +
  geom_col(show.legend = FALSE) + theme_bw()+ scale_fill_manual(values=c("#bcb8b1", "#284b63")) + 
  scale_x_discrete(labels=c(paste0(totals_gr$Var2,"\nn=", totals_gr$total)))+
  geom_text(aes(label = paste0(perc*100,"%")), position = position_stack(vjust = 0.5))+
  labs(x = "", y="Fraction of patients", fill = "me-ctDNA detection status") + 
  theme(plot.title = element_text(hjust = 0.5, face ="bold"), legend.position = "bottom")
  


##### Figure 3H:  Mutations-based ctDNA detection per group ##### 

# BARPLOT 
Barplot_data_HansenData_MUT <- as.data.frame(table(C3_mCRPC_methdata_clinFU$mectDNA_detection, C3_mCRPC_methdata_clinFU$ctDNA_detection))
# ctDNA 
pvalue_HansenData_MUT <- fisher.test(table(C3_mCRPC_methdata_clinFU$mectDNA_detection, C3_mCRPC_methdata_clinFU$ctDNA_detection), workspace = 2e8)$p
totals_HansenData_MUT<- Barplot_data_HansenData_MUT %>%  group_by(Var2) %>%  summarise(total=sum(Freq))

#Percentage
ggplot(Barplot_data_HansenData_MUT %>% group_by(Var2) %>%mutate(perc = round(Freq/sum(Freq),2)), aes(x = Var2, y = perc, fill = Var1, cumulative = TRUE)) +
  geom_col(show.legend = FALSE) + theme_bw()+ scale_fill_manual(values=c("#1F78B4", "#E31A1C")) + 
  geom_text(aes(label = paste0(perc*100,"%")), position = position_stack(vjust = 0.5))+
  labs(x = "", y="Fraction of patients", fill = "me-ctDNA detection status") + 
  theme(plot.title = element_text(hjust = 0.5, face ="bold"), legend.position = "bottom") + scale_x_discrete(labels=c(paste0(totals_HansenData_MUT$Var2,"\nn=", totals_HansenData_MUT$total)))
# eps W270 H225


##### Table S4: Compare sensitivities ##### 

### Compare sensitivity 
table(C3_mCRPC_methdata_clinFU$mectDNA_detection)
C3_mCRPC_detection_table <- table(C3_mCRPC_methdata_clinFU$mectDNA_detection, C3_mCRPC_methdata_clinFU$ctDNA_detection)
fisher.test(C3_mCRPC_detection_table)# fishers exact test, works eventhough there is less than 5 observations in the 2x2 table

library(exact2x2)
mcnemar.exact(C3_mCRPC_detection_table)



#### COHORT 2 SURVIVAL ANALYSES - mCRPC patients only ####### 

# Dichotomize variables  
C3_mCRPC_methdata_clinFU$DMRs_mean_nrpm_GR <- cut(C3_mCRPC_methdata_clinFU$DMRs_mean_nrpm , breaks=c(-Inf,median(C3_mCRPC_methdata_clinFU$DMRs_mean_nrpm ),+Inf), labels=c(0,1)) 
C3_mCRPC_methdata_clinFU$DMRs_mean_nrpm_GR_cohort1cutoff <- cut(C3_mCRPC_methdata_clinFU$DMRs_mean_nrpm , breaks=c(-Inf,1.51,+Inf), labels=c(0,1)) 

C3_mCRPC_methdata_clinFU$PSA_GR <- cut(C3_mCRPC_methdata_clinFU$PSA , breaks=c(-Inf,median(C3_mCRPC_methdata_clinFU$PSA ),+Inf), labels=c(0,1)) 
C3_mCRPC_methdata_clinFU$ctDNA_GR <- cut(C3_mCRPC_methdata_clinFU$ctDNA , breaks=c(-Inf,median(C3_mCRPC_methdata_clinFU$ctDNA ),+Inf), labels=c(0,1)) 

######___ PFS___ ######

#### ___KM plot
#PFS
res.cox_pfs <- coxph(Surv(prog_psa_mo, prog_status)~mectDNA_detection, data=C3_mCRPC_methdata_clinFU)
cox.zph(res.cox_pfs)
summary(res.cox_pfs)


pfs_frag<-survfit(Surv(prog_psa_mo, prog_status)~mectDNA_detection, data=C3_mCRPC_methdata_clinFU)

ggsurvplot(pfs_frag, data=C3_mCRPC_methdata_clinFU,
           combine=TRUE, 
           risk.table=TRUE, 
           conf.int=FALSE, 
           conf.int.style="ribbon",                     # step or ribbon
           censor=TRUE,                                 #TRUE: show censoring (length of FU)
           tables.theme=theme_cleantable(),
           title = "",
           legend.title = "",legend.labs = c("me-ctDNA negative", "me-ctDNA positive"), #legend.labs = c("Low cfMeCaP methylation", "High cfMeCaP methylation"),
           xlab="Progression-free survival (months)",
           ylab="PFS probability",
           palette=c("dodgerblue3","red2", "darkgreen", "purple", "black", "orange", "pink"),                                # other: Set1, jco,
           pval.method = TRUE,
           pval=TRUE, break.x.by = 3,
           pval.coord = c(9,0.9), xlim = c(0,14.5),  
           pval.method.coord = c(9,1), 
           surv.median.line = "hv")  

# Save eps H400, W325 

#OS 
os_frag<-survfit(Surv(os_mo ,os_status)~DMRs_mean_nrpm_GR, data=C3_mCRPC_methdata_clinFU)

ggsurvplot(os_frag, data=C3_mCRPC_methdata_clinFU,
           combine=TRUE, 
           risk.table=TRUE, 
           conf.int=FALSE, 
           conf.int.style="ribbon",                     # step or ribbon
           censor=TRUE,                                 #TRUE: show censoring (length of FU)
           tables.theme=theme_cleantable(),
           title = "",
           legend.title = "",legend.labs = c("Low cfMeCaP48 methylation", "High cfMeCaP48 methylation"), #legend.labs = c("Low cfMeCaP methylation", "High cfMeCaP methylation"),
           xlab="Overall survival (months)",
           ylab="OS probability",
           palette=c("dodgerblue3","red2", "darkgreen", "purple", "black", "orange", "pink"),                                # other: Set1, jco,
           pval.method = TRUE,
           pval=TRUE, break.x.by = 5,
           pval.coord = c(20,0.9), xlim = c(0,30), 
           pval.method.coord = c(20,1), 
           surv.median.line = "hv")  

# intermediate + high ctDNA% patients only
C3_mCRPC_methdata_clinFU_highctDNA <- C3_mCRPC_methdata_clinFU %>% filter(ctDNA > 2)
C3_mCRPC_methdata_clinFU_highctDNA$ctDNA_GRhigh <- cut(C3_mCRPC_methdata_clinFU_highctDNA$ctDNA , breaks=c(-Inf,30,+Inf), labels=c("2-30% ctDNA","30-100% ctDNA")) 


pfs_frag<-survfit(Surv(prog_psa_mo, prog_status)~ctDNA_GRhigh, data=C3_mCRPC_methdata_clinFU_highctDNA)

ggsurvplot(pfs_frag, data=C3_mCRPC_methdata_clinFU_highctDNA,
           combine=TRUE, 
           risk.table=TRUE, 
           conf.int=FALSE, 
           conf.int.style="ribbon",                     # step or ribbon
           censor=TRUE,                                 #TRUE: show censoring (length of FU)
           tables.theme=theme_cleantable(),
           title = "",
           legend.title = "",legend.labs = c("2-30% ctDNA", "30-100% ctDNA"), #legend.labs = c("Low cfMeCaP methylation", "High cfMeCaP methylation"),
           xlab="Progression-free survival (months)",
           ylab="PFS probability",
           palette=c("dodgerblue3","red2", "darkgreen", "purple", "black", "orange", "pink"),                                # other: Set1, jco,
           pval.method = TRUE,
           pval=TRUE, break.x.by = 3,
           pval.coord = c(9,0.9), xlim = c(0,14.5),  
           pval.method.coord = c(9,1), 
           surv.median.line = "hv")  



#### ___Forest plot - univariate cox - PFS ___#### 
Cohort3_SurvData_forest_PFS <- C3_mCRPC_methdata_clinFU %>%
  transmute(prog_psa_mo, prog_status , cfMeCaP = DMRs_mean_nrpm, ctDNA = ctDNA,
            PSA = PSA, 
            cfMeCaP_dichotomized =factor(DMRs_mean_nrpm_GR, labels = c("Low cfMeCaP methylation", "High cfMeCaP methylation")),
            Mets = factor(Metvol),
            ctDNA_dichotomized =factor(ctDNA_GR, labels = c("Low ctDNA%", "High ctDNA%")),
            PSA_dichotomized = factor(PSA_GR, labels = c("PSA<median", "PSA>median")))

# Test assumptions
res.cox_pfs <- coxph(Surv(prog_psa_mo,prog_status )~Mets , data=Cohort3_SurvData_forest_PFS)
cox.zph(res.cox_pfs)
summary(res.cox_pfs)

#### ___Forest plot -Univariate cox - PFS ___#### 

vars_for_table <-c("cfMeCaP","ctDNA","PSA", "Mets")

univ_formulas <- sapply(vars_for_table, function(x) as.formula(paste('Surv(prog_psa_mo, prog_status)~', x)))
univ_models <- lapply(univ_formulas, function(x){coxph(x,data=Cohort3_SurvData_forest_PFS)})
univ_results <- lapply(univ_models,function(x){return(exp(cbind(coef(x),confint(x))))})
print(forest_model(model_list = univ_models,covariates = vars_for_table,merge_models =T))

p.adjust(p = c( 8.34e-05, 0.0021, 0.00207, 0.003343, 0.000863, 0.019082,  0.678302), method="BH")

#### ___Forest plot - Multivariate cox - PFS ___#### 
print(forest_model(coxph(Surv(prog_psa_mo, prog_status)~cfMeCaP + 
                           PSA+ Mets,  data=Cohort3_SurvData_forest_PFS)))

print(forest_model(coxph(Surv(prog_psa_mo, prog_status)~cfMeCaP + ctDNA,  data=Cohort3_SurvData_forest_PFS)))

### __OS 

Cohort3_SurvData_forest_OS <- C3_mCRPC_methdata_clinFU %>%
  transmute(os_mo ,os_status , cfMeCaP = DMRs_mean_nrpm, ctDNA = ctDNA,
            PSA = PSA, 
            cfMeCaP_dichotomized =factor(DMRs_mean_nrpm_GR, labels = c("Low cfMeCaP methylation", "High cfMeCaP methylation")),
            Mets = factor(Metvol),
            ctDNA_dichotomized =factor(ctDNA_GR, labels = c("Low ctDNA%", "High ctDNA%")),
            PSA_dichotomized = factor(PSA_GR, labels = c("Low PSA", "High PSA")))


res.cox_os <- coxph(Surv(os_mo ,os_status )~cfMeCaP + PSA + Mets , data=Cohort3_SurvData_forest_OS)
cox.zph(res.cox_os)
summary(res.cox_os)


#### ___Forest plot -Univariate cox - OS ___#### 

vars_for_table <-c("cfMeCaP_dichotomized","ctDNA_dichotomized","PSA_level", "Metvol_original")
vars_for_table <-c("cfMeCaP","ctDNA","PSA", "Mets")

univ_formulas <- sapply(vars_for_table, function(x) as.formula(paste('Surv(os_mo ,os_status)~', x)))
univ_models <- lapply(univ_formulas, function(x){coxph(x,data=Cohort3_SurvData_forest_OS)})
univ_results <- lapply(univ_models,function(x){return(exp(cbind(coef(x),confint(x))))})
print(forest_model(model_list = univ_models,covariates = vars_for_table,merge_models =T))
# Save eps W700 x H150

# Adj. p 
p.adjust(p = c( 1.45e-05, 0.0016, 0.00262, 0.020359,0.000471, 0.008439,0.551374), method="BH")
p.adjust(p = c( 0.0001, 0.002, 0.034, 0.559, 0.0001), method="BH")


#### ___Forest plot - Multivariate cox - OS ___#### 
print(forest_model(coxph(Surv(os_mo ,os_status)~cfMeCaP  
                           + ctDNA  
                           # + PSA_level  
                           # + Metvol_original
                         , data=Cohort3_SurvData_forest_OS)))
print(forest_model(coxph(Surv(os_mo ,os_status)~cfMeCaP + PSA + Mets, data=Cohort3_SurvData_forest_OS)))


### Mutation-based ctDNA detection 
Cohort3_SurvData_forest_OS$group <- rep("mCRPC", 72)
Cohort3_SurvData_forest_OS$ctDNA_detection[Cohort3_SurvData_forest_OS$ctDNA_fraction > 2] <- "ctDNA Positive (Mutation)"
Cohort3_SurvData_forest_OS$ctDNA_detection[Cohort3_SurvData_forest_OS$ctDNA_fraction < 2] <- "ctDNA Negative (Mutation)"

Barplot_data_gr <- as.data.frame(table(Cohort3_SurvData_forest_OS$ctDNA_detection, Cohort3_SurvData_forest_OS$group))

# ctDNA 
totals_gr<- Barplot_data_gr %>%
  group_by(Var2) %>%
  summarise(total=sum(Freq))

#Percentage
ggplot(Barplot_data_gr %>% group_by(Var2) %>%mutate(perc = round(Freq/sum(Freq),2)), aes(x = Var2, y = perc, fill = Var1, cumulative = TRUE)) +
  geom_col(show.legend = TRUE) + theme_bw()+ scale_fill_manual(values=c("#bcb8b1", "#284b63")) + 
  scale_x_discrete(labels=c(paste0(totals_gr$Var2,"\nn=", totals_gr$total)))+
  geom_text(aes(label = paste0(perc*100,"%")), position = position_stack(vjust = 0.5))+
  labs(x = "", y="Fraction of patients", fill = "me-ctDNA detection status") + 
  theme(plot.title = element_text(hjust = 0.5, face ="bold"), legend.position = "bottom")



## Hansen et al - clin summary
C3_clinSum <- as.data.frame(read_excel("Cohort 3_Sampleinfo_clinical_updated.xlsx", sheet = "Manus_mCRPC_clintable"))

C3_clinSum %>% summarise(PSA_median=median(PSA),
                         PSA_min=min(PSA),
                         PSA_max=max(PSA))




table(C3_clinSum$Age)
sum(table(C3_clinSum$Age)/sum(table(C3_clinSum$Age))*100)

###############################################################################

########## cfMeCaP dynamics during treatment 

#Load dynamics data
Dynamics <- readRDS("PSA_response_data.rds")

######## dynamics groups vs progression within 12 months 

Dynamics$mectDNA_dynamics_newgroups[Dynamics$mectDNA_dynamics %in% c("No clearance", "Partial clearance")] <- "No-clearance" 
Dynamics$mectDNA_dynamics_newgroups[Dynamics$mectDNA_dynamics %in% c("Full clearance", "Undetected")] <- "Clearance" 

# PSA progression 
Dynamics$PSApgr_within1y <- ifelse(test = Dynamics$prog_psa_status == 1 & Dynamics$PFS_months < 12, yes = "Progression", no = "No progression")

Barplot_data_gr <- as.data.frame(table( Dynamics$PSApgr_within1y, Dynamics$mectDNA_dynamics_newgroups))
pvalue_gr <- fisher.test(table(Dynamics$PSApgr_within1y, Dynamics$mectDNA_dynamics_newgroups), workspace = 2e8)$p

totals_gr<- Barplot_data_gr %>%
  group_by(Var2) %>%
  summarise(total=sum(Freq))

#Percentage
ggplot(Barplot_data_gr %>% group_by(Var2) %>%mutate(perc = round(Freq/sum(Freq),2)), aes(x = Var2, y = perc, fill = Var1, cumulative = TRUE)) +
  geom_col(show.legend = FALSE) + theme_bw()+ 
  scale_fill_manual(values=c("#A6CEE3", "#FB9A99")) + 
  #scale_fill_manual(values=c("#1F78B4", "#A6CEE3", "#FB9A99", "#E31A1C")) + 
  scale_x_discrete(labels=c(paste0(totals_gr$Var2,"\nn=", totals_gr$total)))+
  geom_text(aes(label = paste0(perc*100,"%")), position = position_stack(vjust = 0.5))+
  labs(x = "", y="Fraction of patients", fill = "me-ctDNA detection \n status") + 
  theme(plot.title = element_text(hjust = 0.5, face ="bold"), legend.position = "right") 
# save eps 225 x 250


# Radiographic progression 
Dynamics$prog_radio_status_1year <- ifelse(test = (Dynamics$radioPFS_months < 12 & Dynamics$prog_radio_status == 1), 
                                                          yes = "Progression", no = "No progression")

Barplot_data_gr <- as.data.frame(table(Dynamics$prog_radio_status_1year, Dynamics$mectDNA_dynamics_newgroups ))
pvalue_gr <- fisher.test(table(Dynamics$prog_radio_status_1year, Dynamics$mectDNA_dynamics_newgroups), workspace = 2e8)$p
totals_gr<- Barplot_data_gr %>%
  group_by(Var2) %>%
  summarise(total=sum(Freq))

#Percentage
ggplot(Barplot_data_gr %>% group_by(Var2) %>%mutate(perc = round(Freq/sum(Freq),2)), aes(x = Var2, y = perc, fill = Var1, cumulative = TRUE)) +
  geom_col(show.legend = FALSE) + theme_bw()+ 
  scale_fill_manual(values=c("#A6CEE3",  "#FB9A99" )) + 
  #scale_fill_manual(values=c("#1F78B4", "#A6CEE3", "#FB9A99", "#E31A1C")) + 
  scale_x_discrete(labels=c(paste0(totals_gr$Var2,"\nn=", totals_gr$total)))+
  geom_text(aes(label = paste0(perc*100,"%")), position = position_stack(vjust = 0.5))+
  labs(x = "", y="Fraction of patients", fill = "Progression status") + 
  theme(plot.title = element_text(hjust = 0.5, face ="bold"), legend.position = "right") 
# save eps 225 x 250


###### Swimmer plot - figure 4B #######

C2_Swimmer <- readRDS("C2_Swimmerplot_data.rds")
C2_Swimmer$mectDNA_dynamics_newgroups[C2_Swimmer$mectDNA_dynamics.x %in% c("No clearance", "Partial clearance")] <- "No-clearance" 
C2_Swimmer$mectDNA_dynamics_newgroups[C2_Swimmer$mectDNA_dynamics.x %in% c("Full clearance", "Undetected")] <- "Clearance" 

C2_Swimmer$mectDNA_dynamics_newgroups_fac <-  factor(C2_Swimmer$mectDNA_dynamics_newgroups, levels = c("No-clearance", "Clearance"))

#Last clin update
C2_Swimmer$prog_psa_days[C2_Swimmer$PT_NR == 9] <- 126
C2_Swimmer$prog_psa_status[C2_Swimmer$PT_NR == 9] <- 1
C2_Swimmer$PFS_months[C2_Swimmer$PT_NR == 9] <- 4.14

C2_Swimmer$prog_psa_days[C2_Swimmer$PT_NR == 2] <- 80
C2_Swimmer$prog_psa_status[C2_Swimmer$PT_NR == 2] <- 1
C2_Swimmer$PFS_months[C2_Swimmer$PT_NR == 2] <- 2.63



# 
C2_Swimmer %>% 
  ggplot(aes(x=PT_NR.TxFac)) + 
  geom_bar(aes(y=(Enza_tx_max400d_days), fill = factor(tx_ended)), stat = "identity", width = 0.20) + #factor(mectDNA_dynamics)
  scale_y_continuous(expand = c(0,10), limits = c(0,405)) + geom_hline(yintercept = 365, linetype = 2) + 
  coord_flip() + 
  geom_point(aes(y=tt1stscan_days), col = "grey", shape = 3, size = 4, stroke = 2) + 
  geom_point(aes(y=tt2ndscan_days), col = "grey", shape = 3, size = 4, stroke = 2) + 
  geom_point(aes(y=time, shape = MectDNA_detection), col= "black", na.rm = T, size =4, stroke = 1.5) + #col = timepoint_new
  geom_point(aes(y=prog_psa_days), col = "orange", shape = 4, size = 4, stroke = 2) + 
  geom_point(aes(y=prog_radio_days), col = "red", shape = 4, size = 4, stroke = 2) + 
  
  theme_classic() + 
  #scale_fill_manual(labels = c("Undetected", "Full clearance", "Partial clearance", "No clearance"), values=c("#1F78B4", "#A6CEE3", "#FB9A99", "#E31A1C")) + 
  scale_fill_manual(labels = c("Active enzalutamide treatment", "Enzalutamide treatment ended"), values = c( "1" = "#737373", "0" = "#D9D9D9")) +
  
  scale_color_manual(name = "Time point", labels = c("Baseline", "Scan1", "Scan2"),
                     values = c("Baseline" = 'Black',
                                "Scan1" = "#08519C",
                                "Scan2"="#006D2C")) + 
  
  scale_shape_manual(name = "cfMeCaP detection", labels = c( "cfMeCaP_negative", "cfMeCaP_positive"),
                     values = c("cfMeCaP_negative" = 1,
                                "cfMeCaP_positive"=16)) + 
  
  theme(axis.title.x = element_text(face = "bold", size = 14), 
        axis.title.y = element_text(face = "bold", size = 14, vjust = 2)) + 
  facet_grid(rows = vars(mectDNA_dynamics_newgroups_fac), scales = "free_y", space="free_y" ) + 
  labs(title = "cfMeCaP me-ctDNA detection", y= "Time from enzalutamide initiation (days)", x="Patient ID", fill  ="Treatment status")


##################################### Additional data ##########################################

########  Prepare for 450K data analysis ######## 

#Annotate 450K data and extract CpGs of interest 
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
anno <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
anno <- anno[c("chr", "pos", "Name")]
anno_df <- as.data.frame(anno)
anno_df$end <- anno_df$pos + 0
anno_GR <- makeGRangesFromDataFrame(df = anno_df, keep.extra.columns = TRUE, seqnames.field = "chr", start.field = "pos", end.field = "end")

MethSigRegs_GR <- makeGRangesFromDataFrame(df = MethSigRegs, seqnames.field = "chr", start.field = "window_start", end.field = "window_end", keep.extra.columns = T)


hits <- findOverlaps(anno_GR, MethSigRegs_GR,ignore.strand=TRUE)
ovpairs <- Pairs(anno_GR, MethSigRegs_GR, hits=hits)
pint <- pintersect(ovpairs,ignore.strand=TRUE)
mcols(pint) <- data.frame(ovpairs@first, ovpairs@second) #If not running: detach dplyr package. 
pint_df <- as.data.frame(pint)

RegionToCpG <- as.data.frame(pint_df %>% 
                               dplyr::select(region, seqnames.1,start.1, end.1, Name) %>%
                               group_by(seqnames.1, start.1, end.1) %>% 
                               mutate(CpGs = paste(Name, collapse = ";"), .keep="unused") %>%
                               mutate(across(where(is.numeric), mean))%>%
                               distinct(seqnames.1, start.1, .keep_all = TRUE))


RegionToCpG$CpG_counts <- 1+str_count(string = RegionToCpG$CpGs, pattern = ";")

I450K_CpGs_in_signature <- pint_df$Name


##################################################

####  MarmalAid Data ######
Prostate_healthy_450K <- readRDS("FilteredData_Prostate_healthy.rds")
Prostate_healthy_450K_sub <- as.data.frame(t(subset(Prostate_healthy_450K, rownames(Prostate_healthy_450K) %in% I450K_CpGs_in_signature)))
Prostate_healthy_450K_sub$sample_group <- rep("Non-malignant Prostate", length(rownames(Prostate_healthy_450K_sub)))

PC_450K <- readRDS("FilteredData_Prostate_cancer.rds")
PC_450K_sub <- as.data.frame(t(subset(PC_450K, rownames(PC_450K) %in% I450K_CpGs_in_signature)))
PC_450K_sub$sample_group <- rep("Prostate Tumor", length(rownames(PC_450K_sub)))

Blood_450K <- readRDS("FilteredData_Blood.rds")
Blood_450K_sub <- as.data.frame(t(subset(Blood_450K, rownames(Blood_450K) %in% I450K_CpGs_in_signature)))
Blood_450K_sub$sample_group <- rep("Blood", length(rownames(Blood_450K_sub)))


# combine datasets 
PC_Healthy_450K_combined <- rbind(Prostate_healthy_450K_sub, 
                                  PC_450K_sub[colnames(Prostate_healthy_450K_sub)],  
                                  Blood_450K_sub[colnames(Prostate_healthy_450K_sub)]
                                  )


# average probe methylation levels into 300 bp bins. 
CpGs_firstreg <- RegionToCpG$CpGs[RegionToCpG$region == MethSigRegs$region[1]]
CpG_length_firstreg <- length(CpGs_firstreg)
Meth450K_data_firstreg<- PC_Healthy_450K_combined[c(CpGs_firstreg, "sample_group")]

Meth450K_data_firstreg$patient_mean <- apply(X = Meth450K_data_firstreg[c(1:CpG_length_firstreg)],1, mean)

# create df 
M450K_MeanRegMeth_PerPatient <- Meth450K_data_firstreg[2:3]
colnames(M450K_MeanRegMeth_PerPatient)[2] <- MethSigRegs$region[1]


for (i in 2:length(MethSigRegs$region)) {
  
  CpGs <- RegionToCpG$CpGs[RegionToCpG$region == MethSigRegs$region[i]]
  CpG_vect <- unlist(str_split(string = CpGs, pattern = ";"))
  CpG_length <- length(CpG_vect)
  
  gene <- MethSigRegs$newgene[MethSigRegs$region == MethSigRegs$region[i]]
  
  Meth450K_data_sub <- PC_Healthy_450K_combined[c(CpG_vect,"sample_group")]
  
  Meth450K_data_sub$patient_mean <- apply(X = Meth450K_data_sub[c(1:CpG_length)],1, mean)
  
  M450K_MeanRegMeth_PerPatient <- cbind(M450K_MeanRegMeth_PerPatient, reg = Meth450K_data_sub$patient_mean)
  colnames(M450K_MeanRegMeth_PerPatient)[i+1] <- MethSigRegs$region[i]
  
}

# Calculate cfMeCaP score for each patient 
M450K_MeanRegMeth_PerPatient$Mean_Patient_Meth_allregs <- 
  apply(X = M450K_MeanRegMeth_PerPatient[2:length(colnames(M450K_MeanRegMeth_PerPatient))]
        ,1, mean)

my_comparisons_450K <- list(c("Healthy Prostate", "Prostate Tumor"),
                            c( "Blood", "Prostate Tumor"),
                            c("Blood","Healthy Prostate"))


#######  Marmalaid cfMeCap methylation levels 
ggplot(M450K_MeanRegMeth_PerPatient, aes(x = sample_group , y = Mean_Patient_Meth_allregs)) +
  geom_boxplot(aes(fill = sample_group),outlier.shape = NA, show.legend = F, alpha = 1) +
  geom_jitter(position=position_jitter(0.3),  size = 3, alpha=1) +
    scale_fill_manual(values = c("#08519C","#FF7F00", "#A50F15", "#67000D")) +
    theme_classic() +
  labs(title= "" ,y="cfMeCaP48 methylation (nrpm)", x = "") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
        axis.text.x = element_text(size = 12, colour = "black"),
        axis.text.y = element_text(size = 10, colour = "black"),
        axis.title.y = element_text(face = "bold", size = 12, hjust = .5)) + 
  stat_compare_means(comparisons = my_comparisons_450K, method = "wilcox.test", size = 3) + ylim(c(0,1))


# only tissue samples 
M450K_MeanRegMeth_PerPatient %>% filter(sample_group != "Blood") %>% 
ggplot(., aes(x = sample_group , y = Mean_Patient_Meth_allregs)) +
  geom_boxplot(aes(fill = sample_group),outlier.shape = NA, show.legend = F, alpha = 1) +
  geom_jitter(position=position_jitter(0.3),  size = 3, alpha=1) +
  scale_fill_manual(values = c("#08519C", "#A50F15")) +
  theme_classic() +
  labs(title= "" ,y="cfMeCaP48 methylation (nrpm)", x = "") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
        axis.text.x = element_text(size = 12, colour = "black"),
        axis.text.y = element_text(size = 10, colour = "black"),
        axis.title.y = element_text(face = "bold", size = 12, hjust = .5)) + 
  stat_compare_means(method = "wilcox.test", size = 3) + ylim(c(0,1))


M450K_MeanRegMeth_PerPatient %>% filter(sample_group != "Blood") %>% 
  ggplot(aes(x = sample_group, y = Mean_Patient_Meth_allregs, fill = sample_group)) +
  geom_violin(alpha = 1, show.legend = F, scale = "width") +
  geom_boxplot(width = 0.15, outlier.shape = NA, show.legend = F, 
               fill = "white", alpha = 1) +
  scale_fill_manual(values = c("#08519C", "#FED976")) +
  scale_color_manual(values = c("#08519C", "#FED976")) +
  labs(title= "" ,y="cfMeCaP methylation (beta)", x = "") +
  theme_classic() + coord_cartesian(ylim = c(0,1)) + 
  theme(axis.text.x = element_text(size = 12, colour = "black"),
        axis.text.y = element_text(size = 10, colour = "black"),
        axis.title.y = element_text(face = "bold", size = 12, hjust = 0.5)) +
  stat_compare_means(method = "wilcox.test", size = 3)



####### EPIC data #####

## annotate EPIC methylation data 
library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)

anno_EIPC <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
anno_EIPC <- anno_EIPC[c("chr", "pos", "Name")]
anno_EPIC_df <- as.data.frame(anno_EIPC)
anno_EPIC_df$end <- anno_EPIC_df$pos + 0
 
## Load GSE152026 data 
EUGEI_processed_signals_ChAMP_normalized <- readRDS("GSE152026_EUGEI_processed_signals_ChAMP_normalized.rds")

# Annotate data
EUGEI_annotated <- cbind(EUGEI_processed_signals_ChAMP_normalized, 
                         anno_EPIC_df[rownames(EUGEI_processed_signals_ChAMP_normalized),])

# Create GR object
anno_EPIC_GR <- makeGRangesFromDataFrame(df = anno_EPIC_df, keep.extra.columns = TRUE, seqnames.field = "chr", start.field = "pos", end.field = "end")


# Identify cfMeCaP overlapping regions
hits_EPIC <- findOverlaps(anno_EPIC_GR, MethSigRegs_GR,ignore.strand=TRUE)
ovpairs_EPIC <- Pairs(anno_EPIC_GR, MethSigRegs_GR, hits=hits_EPIC)
pint_EPIC <- pintersect(ovpairs_EPIC,ignore.strand=TRUE)
mcols(pint_EPIC) <- data.frame(ovpairs_EPIC@first, ovpairs_EPIC@second) #If not running: detach dplyr package. 
pint_EPIC_df <- as.data.frame(pint_EPIC)


RegionToCpG_EPIC <- as.data.frame(pint_EPIC_df %>% 
                                    dplyr::select(region, seqnames.1,start.1, end.1, Name) %>%
                                    group_by(seqnames.1, start.1, end.1) %>% 
                                    mutate(CpGs = paste(Name, collapse = ";"), .keep="unused") %>%
                                    mutate(across(where(is.numeric), mean))%>%
                                    distinct(seqnames.1, start.1, .keep_all = TRUE))

EPIC_CpGs_in_signature <- pint_EPIC_df$Name


# Two cfMeCaP regions does not overlap the EPIC data. 
# Remove these regions to compare data 
MethSigRegs_EPICcomp <- subset(MethSigRegs, region %in% RegionToCpG_EPIC$region)

## Extract signature CpG sites from EPIC data 
EUGEI_sub <- as.data.frame(t(subset(EUGEI_processed_signals_ChAMP_normalized, rownames(EUGEI_processed_signals_ChAMP_normalized) %in% EPIC_CpGs_in_signature)))

EUGEI_CpGs <- colnames(EUGEI_sub)
EPIC_CpGs_length <- length(EUGEI_CpGs)

# Calculate mean methylation per patient 
EPIC_CpGs_firstreg <- RegionToCpG_EPIC$CpGs[RegionToCpG_EPIC$region == MethSigRegs_EPICcomp$region[1]]
EPIC_CpG_length_firstreg <- length(EPIC_CpGs_firstreg)
EUGEI_sub_firstreg<- EUGEI_sub[c(EPIC_CpGs_firstreg)]

EUGEI_sub_firstreg$patient_mean <- apply(X = EUGEI_sub_firstreg[c(1:EPIC_CpG_length_firstreg)],1, mean)

# create df 
EUGEI_MeanRegMeth_PerPatient <- EUGEI_sub_firstreg[2]
colnames(EUGEI_MeanRegMeth_PerPatient)[1] <- MethSigRegs_EPICcomp$region[1]

for (i in 2:length(MethSigRegs_EPICcomp$region)) {
  
  EPIC_CpGs <- RegionToCpG_EPIC$CpGs[RegionToCpG_EPIC$region == MethSigRegs_EPICcomp$region[i]]
  EPIC_CpG_vect <- unlist(str_split(string = EPIC_CpGs, pattern = ";"))
  EPIC_CpG_vect <- EPIC_CpG_vect[EPIC_CpG_vect %in% EUGEI_CpGs]
  
  EPIC_CpG_length <- length(EPIC_CpG_vect)
  
  gene <- MethSigRegs_EPICcomp$newgene[MethSigRegs_EPICcomp$region == MethSigRegs_EPICcomp$region[i]]
  
  EUGEI_data_sub <- EUGEI_sub[c(EPIC_CpG_vect)]
  
  EUGEI_data_sub$patient_mean <- apply(X = EUGEI_data_sub[c(1:EPIC_CpG_length)],1, mean)
  
  EUGEI_MeanRegMeth_PerPatient <- cbind(EUGEI_MeanRegMeth_PerPatient, reg = EUGEI_data_sub$patient_mean)
  colnames(EUGEI_MeanRegMeth_PerPatient)[i] <- MethSigRegs_EPICcomp$region[i]
  
}


EUGEI_MeanRegMeth_PerPatient


#### Calculate cfMeCaP48 score for each patient 

EUGEI_MeanRegMeth_PerPatient$Mean_Patient_Meth_allregs <- 
  apply(X = EUGEI_MeanRegMeth_PerPatient[2:length(colnames(EUGEI_MeanRegMeth_PerPatient))]
        ,1, mean)

EUGEI_MeanRegMeth_PerPatient$group <- rep("Blood", length(rownames(EUGEI_MeanRegMeth_PerPatient)))

#### Load sample info 
#setwd("E:/molpros/faststorage/BACKUP/Karina/Methylation_Data_450K/GEO data/EPIC_methylation data/GSE152026 (Blood from psychotic patients + healthy controls, EPIC)")
EPIC_sampleinfo <- as.data.frame(read_xlsx(path = "GSE152026_sampleinfo.xlsx", sheet = "Selected"))


EPIC_sampleinfo$SampleID <- str_remove(string = EPIC_sampleinfo$Sample_title, pattern = " Sample*.")
EPIC_ctrls <- subset(EPIC_sampleinfo,Phenotype == "Control") 
EPIC_ctrls$Age_nu <- as.numeric(EPIC_ctrls$Age)
length(EPIC_ctrls$Sample_title)# 245 controls
summary(EPIC_ctrls$Age_nu)

summary(EUGEI_MeanRegMeth_PerPatient$Mean_Patient_Meth_allregs)

#boxplot
ggplot(data = EUGEI_MeanRegMeth_PerPatient, aes(x=group, y = Mean_Patient_Meth_allregs)) +
  geom_boxplot(aoutlier.shape = NA, show.legend = F, alpha = 1, color ="#08519C") +
  geom_jitter(position=position_jitter(0.5),  size = 4, alpha=1) +
  theme_classic() + ylim(0,1) +
  labs(       y="Mean signature methylation (beta)", x = "") +
  theme(axis.text.x = element_text(size = 10, colour = "black"),
        axis.text.y = element_text(size = 10, colour = "black"),
        axis.title.y = element_text(face = "bold", size = 10))

#violin
  ggplot(data = EUGEI_MeanRegMeth_PerPatient, aes(x=group, y = Mean_Patient_Meth_allregs, fill = group)) +
  geom_violin(alpha = 1, show.legend = F, scale = "width") +
  geom_boxplot(width = 0.15, outlier.shape = NA, show.legend = F, 
               fill = "white", alpha = 1) +
  scale_fill_manual(values = c("#08519C")) +
  scale_color_manual(values = c("#08519C")) +
  labs(title= "" ,y="cfMeCaP methylation (beta)", x = "") +
  theme_classic() + coord_cartesian(ylim = c(0,1)) + 
  theme(axis.text.x = element_text(size = 12, colour = "black"),
        axis.text.y = element_text(size = 10, colour = "black"),
        axis.title.y = element_text(face = "bold", size = 12, hjust = 0.5)) +
  stat_compare_means(method = "wilcox.test", size = 3)



