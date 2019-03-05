setwd("~/Desktop/Oxford Third Year Work/FHS Project/R/Data Analysis/Analysis of Cleaned Up Data")
library(ggplot2)
library(cowplot)
library(DescTools)
library(EnvStats)
library(tidyr)
library(rlang)
library(stringi)


Ct_values <- read.csv("Ct Values Cleaned.csv")


#Boxplot of Ct Values - to look at raw spread of each primer

Ct_values_boxplot <- ggplot(Ct_values, aes(x=Primer, y=Mean.Ct, color=Primer)) +
  geom_boxplot() +
  ggtitle("Mean Ct Values") +
  scale_x_discrete(labels=c("miR122", "miR128", "miR223", "miR24", "miR27", "miR29", "miR320", "miR378", "miR4497"))


#long to wide format

Ct_values_wide <- reshape(Ct_values, idvar = "Participant.ID", timevar = "Primer", direction = "wide")
 
rownames(Ct_values_wide) <- Ct_values_wide[ ,1]
Ct_values_wide <- Ct_values_wide[ ,-1]


#calculating means of different normaliser groups

miR_normaliser_means <- apply(Ct_values_wide[ ,c(1,3)], 1, geoMean)

total_normaliser_means <- apply(Ct_values_wide[ ,5:6], 1, geoMean)

normfinder_means <- apply(Ct_values_wide[ ,c(2,6)], 1, geoMean)

#df of just target primers

target_Cts <- Ct_values_wide[ ,7:9]
colnames(target_Cts) <- c("miR122", "miR4497", "miR223")


#making dfs of just pre and just post Cts

pre_samples <- grep("A", rownames(Ct_values_wide))

pre_cts <- Ct_values_wide[pre_samples, ]

post_samples <- grep("B", rownames(Ct_values_wide))

post_cts <- Ct_values_wide[post_samples, ]


delta_cts <- post_cts - pre_cts


#boxplot of delta Cts

?gather

long_delta_cts <- gather(delta_cts, Primer, Delta.Ct, 1:9)

Delta_Ct_boxplot <- ggplot(long_delta_cts, aes(x=Primer, y=Delta.Ct, color=Primer)) +
  geom_boxplot() +
  scale_x_discrete(labels=c("miR122", "miR128", "miR223", "miR24", "miR27", "miR29", "miR320", "miR378", "miR4497"))


#CIRCULAR ANALYSIS

circular_analysis_delta_cts <- delta_cts[ ,1:3]

c_a_320_128 <- circular_analysis_delta_cts[ ,2]-circular_analysis_delta_cts[ ,1]
c_a_378_320 <- circular_analysis_delta_cts[ ,3]-circular_analysis_delta_cts[ ,2]
c_a_128_378 <- circular_analysis_delta_cts[ ,1]-circular_analysis_delta_cts[ ,3]

c_a_320_128 <- data.frame(c_a_320_128)
c_a_378_320 <- data.frame(c_a_378_320)
c_a_128_378 <- data.frame(c_a_128_378)

circular_analysis_DDCT <- cbind(c_a_320_128, c_a_378_320, c_a_128_378)

colnames(circular_analysis_DDCT) <- c("miR320 vs miR128", "miR378 vs miR320", "miR128 vs miR378")

c_a_wilcox <- apply(circular_analysis_DDCT, 2, wilcox.test, conf.int=TRUE)
 

#taking miRs normaliser means away from target_cts

miR_normalised_cts <- target_Cts - miR_normaliser_means


pre_miR_normalised_Cts <- miR_normalised_cts[pre_samples, ]
post_miR_normalised_Cts <- miR_normalised_cts[post_samples, ]

miR_DDCT <- post_miR_normalised_Cts - pre_miR_normalised_Cts
colnames(miR_DDCT) <- c("miR-122", "miR-4497", "miR-223")

long_miR_DDCT <- gather(miR_DDCT, Target_Gene, Delta.delta.Ct, 1:3)

miR_DDCT_boxplot <- ggplot(long_miR_DDCT, aes(x=Target_Gene, y=-Delta.delta.Ct, color=Target_Gene)) +
  geom_boxplot()+
  ylab("Log2 Fold Change")+
  xlab("Target Gene")+
  scale_x_discrete(labels=c("miR-122", "miR-223", "miR-4497"))+
  geom_hline(yintercept=0, linetype="dashed")


miR_wilcox <- apply(miR_DDCT, 2, wilcox.test, conf.int=TRUE)



#using normfinder normalisers

normfinder_normalised_cts <- target_Cts - normfinder_means


pre_normfinder_normalised_Cts <- normfinder_normalised_cts[pre_samples, ]
post_normfinder_normalised_Cts <- normfinder_normalised_cts[post_samples, ]

normfinder_DDCT <- post_normfinder_normalised_Cts - pre_normfinder_normalised_Cts

long_normfinder_DDCT <- gather(normfinder_DDCT, Primer, Delta.delta.Ct, 1:3)

normfinder_DDCT_boxplot <- ggplot(long_normfinder_DDCT, aes(x=Primer, y=-Delta.delta.Ct, color=Primer)) +
  geom_boxplot() +
  ggtitle("Log2 Fold Change: Normfinder Normalisers")+
  ylab("Log2 Fold Change")+
  scale_x_discrete(labels=c("miR122", "miR223", "miR4497"))


normfinder_wilcox <- apply(normfinder_DDCT, 2, wilcox.test, conf.int=TRUE)


#using total normalisers

total_normalised_cts <- target_Cts - total_normaliser_means


pre_total_normalised_Cts <- total_normalised_cts[pre_samples, ]
post_total_normalised_Cts <- total_normalised_cts[post_samples, ]

total_DDCT <- post_total_normalised_Cts - pre_total_normalised_Cts

long_total_DDCT <- gather(total_DDCT, Primer, Delta.delta.Ct, 1:3)

total_DDCT_boxplot <- ggplot(long_total_DDCT, aes(x=Primer, y=-Delta.delta.Ct, color=Primer)) +
  geom_boxplot()+
  ggtitle("Log2 Fold Change: Total Normalisers")+
  ylab("Log2 Fold Change")+
  scale_x_discrete(labels=c("miR122", "miR223", "miR4497"))


total_wilcox <- apply(total_DDCT, 2, wilcox.test, conf.int=TRUE)


#Combined log2 fold change boxplots

Log2_FC <- plot_grid(miR_DDCT_boxplot, total_DDCT_boxplot, normfinder_DDCT_boxplot)



#Routine vs MenB�

routine_MenB <- c("Routine", "Routine", "Routine", "Routine", "Routine", "Routine", "MenB", "MenB", "MenB", "Routine", "MenB", "Routine", "Routine", "MenB", "MenB", "Routine", "Routine", "MenB", "MenB", "MenB")

R_M_miR_DDCT <- cbind(miR_DDCT, routine_MenB)
R_M_normfinder_DDCT <- cbind(normfinder_DDCT, routine_MenB)
R_M_total_DDCT <- cbind(total_DDCT, routine_MenB)


#routine vs MenB miR normalisers

long_R_M_miR_DDCT <- gather(R_M_miR_DDCT, Primer, DDCT, "miR-122":"miR-223")
long_R_M_normfinder_DDCT <- gather(R_M_normfinder_DDCT, Primer, DDCT, miR122:miR223)
long_R_M_total_DDCT <- gather(R_M_total_DDCT, Primer, DDCT, miR122:miR223)

colnames(long_R_M_miR_DDCT) <- c("Group", "Target_Gene", "DDCT")

R_M_miR_DDCT_boxplot <- ggplot(long_R_M_miR_DDCT, aes(x=Target_Gene, y=-DDCT, color=Group))+
  geom_boxplot(aes(group=interaction(Group, Target_Gene)))+
  ylab("Log2 Fold Change")+
  geom_hline(yintercept=0, linetype="dashed")

R_M_total_DDCT_boxplot <- ggplot(long_R_M_total_DDCT, aes(x=Primer, y=-DDCT, color=routine_MenB))+
  geom_boxplot(aes(group=interaction(routine_MenB, Primer)))+
  ggtitle("Log2 Fold Change: Routine vs MenB, total normalisers")+
  ylab("Log2 Fold Change")

R_M_normfinder_DDCT_boxplot <- ggplot(long_R_M_normfinder_DDCT, aes(x=Primer, y=-DDCT, color=routine_MenB))+
  geom_boxplot(aes(group=interaction(routine_MenB, Primer)))+
  ggtitle("Log2 Fold Change: Routine vs MenB, normfinder normalisers")+
  ylab("Log2 Fold Change")

Routine_MenB_Log2_FC <- plot_grid(R_M_miR_DDCT_boxplot, R_M_total_DDCT_boxplot, R_M_normfinder_DDCT_boxplot)


#Wilcoxon signed rank test of routine vs MenB

routine_samples <- grep("Routine", R_M_miR_DDCT[ , 4])
MenB_samples <- grep("MenB", R_M_miR_DDCT[ , 4])

R_M_miR_routine <- R_M_miR_DDCT[routine_samples, ]
R_M_miR_MenB <- R_M_miR_DDCT[MenB_samples, ]

R_M_total_routine <- R_M_total_DDCT[routine_samples, ]
R_M_total_MenB <- R_M_total_DDCT[MenB_samples, ]

R_M_normfinder_routine <- R_M_normfinder_DDCT[routine_samples, ]
R_M_normfinder_MenB <- R_M_normfinder_DDCT[MenB_samples, ]


routine_miR_wilcox <- apply(R_M_miR_routine[ ,1:3], 2, wilcox.test, conf.int=TRUE)
MenB_miR_wilcox <- apply(R_M_miR_MenB[ ,1:3], 2, wilcox.test, conf.int=TRUE)

routine_total_wilcox <- apply(R_M_total_routine[ ,1:3], 2, wilcox.test, conf.int=TRUE)
MenB_total_wilcox <- apply(R_M_total_MenB[ ,1:3], 2, wilcox.test, conf.int=TRUE)

routine_normfinder_wilcox <- apply(R_M_normfinder_routine[ ,1:3], 2, wilcox.test, conf.int=TRUE)
MenB_normfinder_wilcox <- apply(R_M_normfinder_MenB[ ,1:3], 2, wilcox.test, conf.int=TRUE)


#correlating PCR data with sequencing data

sequencing_data <- read.csv("Sequencing Data.csv")

primers <- c("hsa-miR-128-3p","hsa-miR-320a","hsa-miR-378a-3p","hsa-miR-24-3p","hsa-miR-27a-3p","hsa-miR-29a-3p","hsa-miR-122-5p","hsa-miR-4497","hsa-miR-223-5p")

primers

rownames(sequencing_data) <- sequencing_data[ ,1] 
sequencing_data <- sequencing_data[ ,-1]

#how to find a list within a list - but names must match exactly
#taking mimat number out of row names so that row name matches 'hsa-miR-320a' etc exactly

sequencing_data2 <- do.call(rbind, strsplit(rownames(sequencing_data), ' '))  

rownames(sequencing_data) <- sequencing_data2[,1]

true <- rownames(sequencing_data)%in%primers #gives vectors of true and false

sequencing_data <- sequencing_data[true, ] #keeps rows that are true

log2_sequencing_data <- log2(sequencing_data)

log2_sequencing_data <- log2_sequencing_data[ ,c(-1, -34)]

Participant_ID <- c("97A", "25B", "91A", "23A", "103B", "052B", "149B", "103A", "152B", "88B", "052A", "61B", "23B", "123A", "88A", "97B", "133A", "33B", "19B", "123B", "8087B", "143A", "152A", "143B", "33A", "91B", "8892B", "25A", "149A", "27A", "112B", "19A", "112A", "133B", "8123A", "27B", "8087A", "8892A", "8123B", "61A")

colnames(log2_sequencing_data)[1:40] <- Participant_ID

#joining each of the reference groups into separate data frames in order to average them


miR_references <- c("hsa-miR-128-3p","hsa-miR-378a-3p")
true_miRs <- rownames(log2_sequencing_data)%in%miR_references
sequencing_miR_references_df <- log2_sequencing_data[true_miRs, ]

total_references <- c("hsa-miR-24-3p","hsa-miR-27a-3p","hsa-miR-29a-3p")
true_total <- rownames(log2_sequencing_data)%in%total_references
sequencing_total_references_df <- log2_sequencing_data[true_total, ]

normfinder_references <- c("hsa-miR-320a", "hsa-miR-29a-3p")
true_normfinder <- rownames(log2_sequencing_data)%in%normfinder_references
sequencing_normfinder_references_df <- log2_sequencing_data[true_normfinder, ]

#Calculating geometric means

miR_reference_geomean <- apply(sequencing_miR_references_df, 2, Gmean)   #applies function 'Gmean' to 2nd dimension (columns)
sequencing_miR_references_df <- rbind(sequencing_miR_references_df, miR_reference_geomean)
rownames(sequencing_miR_references_df)[3] <- "Geometric Mean"

total_reference_geomean <- apply(sequencing_total_references_df, 2, Gmean)   #applies function 'Gmean' to 2nd dimension (columns)
sequencing_total_references_df <- rbind(sequencing_total_references_df, total_reference_geomean)
rownames(sequencing_total_references_df)[4] <- "Geometric Mean"

normfinder_reference_geomean <- apply(sequencing_normfinder_references_df, 2, Gmean)   #applies function 'Gmean' to 2nd dimension (columns)
sequencing_normfinder_references_df <- rbind(sequencing_normfinder_references_df, normfinder_reference_geomean)
rownames(sequencing_normfinder_references_df)[3] <- "Geometric Mean"

#making a data frame of target miRNAs

target_miRNAs <- c("hsa-miR-122-5p","hsa-miR-4497","hsa-miR-223-5p")
true_targets <- rownames(log2_sequencing_data)%in%target_miRNAs
sequencing_target_miRNAs <- log2_sequencing_data[target_miRNAs, ]

sequencing_DeltaCT_miR122_miRs<-sequencing_target_miRNAs[1, ] - sequencing_miR_references_df[3, ]
sequencing_DeltaCT_miR122_total <- sequencing_target_miRNAs[1, ] - sequencing_total_references_df[4, ]
sequencing_DeltaCT_miR122_normfinder <- sequencing_target_miRNAs[1, ] - sequencing_normfinder_references_df[3, ]


sequencing_DeltaCT_miR4497_miRs <- sequencing_target_miRNAs[2, ] - sequencing_miR_references_df[4, ]
sequencing_DeltaCT_miR4497_total <- sequencing_target_miRNAs[2, ] - sequencing_total_references_df[4, ]
sequencing_DeltaCT_miR4497_normfinder <- sequencing_target_miRNAs[1, ] - sequencing_normfinder_references_df[3, ]


sequencing_DeltaCT_miR223_miRs <- sequencing_target_miRNAs[3, ] - sequencing_miR_references_df[4, ]
sequencing_DeltaCT_miR223_total <- sequencing_target_miRNAs[3, ] - sequencing_total_references_df[4, ]
sequencing_DeltaCT_miR223_normfinder <- sequencing_target_miRNAs[1, ] - sequencing_normfinder_references_df[3, ]

miRs_count_diff <- rbind(sequencing_DeltaCT_miR122_miRs, sequencing_DeltaCT_miR4497_miRs, sequencing_DeltaCT_miR223_miRs)
miRs_count_diff <- t(miRs_count_diff)
colnames(miRs_count_diff) <- c("miR122", "miR4497", "miR223")

total_count_diff <- rbind(sequencing_DeltaCT_miR122_total, sequencing_DeltaCT_miR4497_total, sequencing_DeltaCT_miR223_total)
total_count_diff <- t(total_count_diff)
colnames(total_count_diff) <- c("miR122", "miR4497", "miR223")

normfinder_count_diff <- rbind(sequencing_DeltaCT_miR122_normfinder, sequencing_DeltaCT_miR4497_normfinder, sequencing_DeltaCT_miR223_normfinder)
normfinder_count_diff <- t(normfinder_count_diff)
colnames(normfinder_count_diff) <- c("miR122", "miR4497", "miR223")


count_pre_samples <- grep("A", rownames(miRs_count_diff))  
count_post_samples <- grep("B", rownames(miRs_count_diff))

miRs_count_pre_normalised <- miRs_count_diff[count_pre_samples, ]
miRs_count_post_normalised <- miRs_count_diff[count_post_samples, ]

total_count_pre_normalised <- total_count_diff[count_pre_samples, ]
total_count_post_normalised <- total_count_diff[count_post_samples, ]

normfinder_count_pre_normalised <- normfinder_count_diff[count_pre_samples, ]
normfinder_count_post_normalised <- normfinder_count_diff[count_post_samples, ]


row_order_pre <- order(rownames(miRs_count_pre_normalised))
miRs_count_pre_normalised <- miRs_count_pre_normalised[row_order_pre, ]
row_order_post <- order(rownames(miRs_count_post_normalised))
miRs_count_post_normalised <- miRs_count_post_normalised[row_order_post, ]

row_order_pre <- order(rownames(total_count_pre_normalised))
total_count_pre_normalised <- total_count_pre_normalised[row_order_pre, ]
row_order_post <- order(rownames(total_count_post_normalised))
total_count_post_normalised <- total_count_post_normalised[row_order_post, ]

row_order_pre <- order(rownames(normfinder_count_pre_normalised))
normfinder_count_pre_normalised <- normfinder_count_pre_normalised[row_order_pre, ]
row_order_post <- order(rownames(normfinder_count_post_normalised))
normfinder_count_post_normalised <- normfinder_count_post_normalised[row_order_post, ]




check_1 <- stri_sub(rownames(miRs_count_pre_normalised), 1, -2)
check_2 <- stri_sub(rownames(miRs_count_post_normalised), 1, -2)


table(check_1 == check_2)


sequencing_log2_FC_miRs <- miRs_count_post_normalised - miRs_count_pre_normalised
sequencing_log2_FC_total <- total_count_post_normalised - total_count_pre_normalised
sequencing_log2_FC_normfinder <- normfinder_count_post_normalised - normfinder_count_pre_normalised


DDCT_order <- order(rownames(miR_DDCT))
miR_DDCT <- miR_DDCT[DDCT_order, ]

DDCT_order <- order(rownames(total_DDCT))
total_DDCT <- total_DDCT[DDCT_order, ]

DDCT_order <- order(rownames(normfinder_DDCT))
normfinder_DDCT <- normfinder_DDCT[DDCT_order, ]


#scatterplots of miR122: miR, total and normfinder normalisers


miR122_miR_log2FC_correlation <- ggplot(mapping=aes(x=sequencing_log2_FC_miRs[ ,1], y=-miR_DDCT[ ,1]))+
  geom_point(colour="red", shape=18, size=2)+
  xlab("Sequencing Log2 Fold Change: miR-122")+
  ylab("PCR Log2 Fold Change: miR-122")+
  geom_smooth(method=lm, colour="indianred1", fill="grey80")

miR122_total_log2FC_correlation <- ggplot(mapping=aes(x=sequencing_log2_FC_total[ ,1], y=-total_DDCT[ ,1]))+
  geom_point(colour="white")+
  xlab("Sequencing log2 Fold Change")+
  ylab("PCR log2 Fold Change")+
  ggtitle("miR122: total normalisers")+ geom_text(label=rownames(miR_DDCT))+
  geom_smooth(method=lm)

miR122_normfinder_log2FC_correlation <- ggplot(mapping=aes(x=sequencing_log2_FC_normfinder[ ,1], y=-normfinder_DDCT[ ,1]))+
  geom_point(colour="white")+
  xlab("Sequencing log2 Fold Change")+
  ylab("PCR log2 Fold Change")+
  ggtitle("miR122: normfinder normalisers")+
  geom_text(label=rownames(miR_DDCT))+
  geom_smooth(method=lm)

miR122_log2FC_correlation <- plot_grid(miR122_miR_log2FC_correlation, miR122_total_log2FC_correlation, miR122_normfinder_log2FC_correlation)


#scatterplots of miR4497: miR, total and normfinder normalisers

miR4497_miR_log2FC_correlation <- ggplot(mapping=aes(x=sequencing_log2_FC_miRs[ ,2], y=-miR_DDCT[ ,2]))+
  geom_point(colour="blue1", shape=18, size=2)+
  xlab("Sequencing log2 Fold Change: miR-4497")+
  ylab("PCR Log2 Fold Change: miR-4497")+
  geom_smooth(method=lm, colour="royalblue2", fill="grey80")

miR4497_total_log2FC_correlation <- ggplot(mapping=aes(x=sequencing_log2_FC_total[ ,2], y=-total_DDCT[ ,2]))+
  geom_point(colour="white")+
  xlab("Sequencing log2 Fold Change")+
  ylab("PCR log2 Fold Change")+
  ggtitle("miR4497: total normalisers")+ geom_text(label=rownames(miR_DDCT))+
  geom_smooth(method=lm)

miR4497_normfinder_log2FC_correlation <- ggplot(mapping=aes(x=sequencing_log2_FC_normfinder[ ,2], y=-normfinder_DDCT[ ,2]))+
  geom_point(colour="white")+
  xlab("Sequencing log2 Fold Change")+
  ylab("PCR log2 Fold Change")+
  ggtitle("miR4497: normfinder normalisers")+
  geom_text(label=rownames(miR_DDCT))+
  geom_smooth(method=lm)

miR4497_log2FC_correlation <- plot_grid(miR4497_miR_log2FC_correlation, miR4497_total_log2FC_correlation, miR4497_normfinder_log2FC_correlation)


#scatterplots of miR223: miR, total and normfinder normalisers

miR223_miR_log2FC_correlation <- ggplot(mapping=aes(x=sequencing_log2_FC_miRs[ ,3], y=-miR_DDCT[ ,3]))+
  geom_point(colour="green4", shape=18, size=2)+
  xlab("Sequencing log2 Fold Change: miR-223")+
  ylab("PCR log2 Fold Change: miR-223")+
  geom_smooth(method=lm, colour="limegreen", fill="grey80")

miR223_total_log2FC_correlation <- ggplot(mapping=aes(x=sequencing_log2_FC_total[ ,3], y=-total_DDCT[ ,3]))+
  geom_point(colour="white")+
  xlab("Sequencing log2 Fold Change")+
  ylab("PCR log2 Fold Change")+
  ggtitle("miR223: total normalisers")+ geom_text(label=rownames(miR_DDCT))+
  geom_smooth(method=lm)

miR223_normfinder_log2FC_correlation <- ggplot(mapping=aes(x=sequencing_log2_FC_normfinder[ ,3], y=-normfinder_DDCT[ ,3]))+
  geom_point(colour="white")+
  xlab("Sequencing log2 Fold Change")+
  ylab("PCR log2 Fold Change")+
  ggtitle("miR223: normfinder normalisers")+
  geom_text(label=rownames(miR_DDCT))+
  geom_smooth(method=lm)

miR223_log2FC_correlation <- plot_grid(miR223_miR_log2FC_correlation, miR223_total_log2FC_correlation, miR223_normfinder_log2FC_correlation)


miR_normaliser_log2FC_correlation <- plot_grid(miR122_miR_log2FC_correlation, miR4497_miR_log2FC_correlation, miR223_miR_log2FC_correlation)


#calculating spearmans rank correlation coefficient for miR-normalised

miR122_miR_corr <- cor.test(x=sequencing_log2_FC_miRs[ ,1], y=-miR_DDCT[ ,1], method = 'spearman')

miR4497_miR_corr <- cor.test(x=sequencing_log2_FC_miRs[ ,2], y=-miR_DDCT[ ,2], method = 'spearman')

miR223_miR_corr <- cor.test(x=sequencing_log2_FC_miRs[ ,3], y=-miR_DDCT[ ,3], method = 'spearman')


#correlation of sequencing vs routine and MenB

sequencing_log2_FC_miRs <- as.data.frame(sequencing_log2_FC_miRs)

count_routine_MenB <- c("Routine", "MenB", "MenB", "Routine", "MenB", "Routine", "Routine", "Routine", "MenB", "MenB", "Routine", "MenB", "Routine", "Routine", "MenB", "Routine", "Routine", "MenB", "MenB", "Routine")

sequencing_log2_FC_miRs <- cbind(sequencing_log2_FC_miRs, count_routine_MenB) 

count_routine_samples <- grep("Routine", sequencing_log2_FC_miRs[ , 4])
count_MenB_samples <- grep("MenB", sequencing_log2_FC_miRs[ , 4])

sequencing_log2_FC_miRs_routine <- sequencing_log2_FC_miRs[count_routine_samples, ]
sequencing_log2_FC_miRs_MenB <- sequencing_log2_FC_miRs[count_MenB_samples, ]

sequencing_log2_FC_miRs_routine <- sequencing_log2_FC_miRs_routine[rownames(R_M_miR_routine),,drop=FALSE]
sequencing_log2_FC_miRs_MenB <- sequencing_log2_FC_miRs_MenB[rownames(R_M_miR_MenB),,drop=FALSE]

sequencing_log2_FC_miRs_routine <- data.frame(sequencing_log2_FC_miRs_routine)
sequencing_log2_FC_miRs_MenB <- data.frame(sequencing_log2_FC_miRs_MenB)
R_M_miR_MenB <- data.frame(R_M_miR_MenB)
R_M_miR_routine <- data.frame(R_M_miR_routine)

#routine correlations

miR122_miR_log2FC_correlation_routine <- ggplot(mapping=aes(x=sequencing_log2_FC_miRs_routine[ ,1], y=-R_M_miR_routine[ ,1]))+
  geom_point(colour="white")+
  xlab("Sequencing log2 Fold Change")+
  ylab("PCR log2 Fold Change")+
  ggtitle("miR122: miRs normalisers, routine")+
  geom_text(label=rownames(R_M_miR_routine))+
  geom_smooth(method=lm)

miR4497_miR_log2FC_correlation_routine <- ggplot(mapping=aes(x=sequencing_log2_FC_miRs_routine[ ,2], y=-R_M_miR_routine[ ,2]))+
  geom_point(colour="white")+
  xlab("Sequencing log2 Fold Change")+
  ylab("PCR log2 Fold Change")+
  ggtitle("miR4497: miRs normalisers, routine")+
  geom_text(label=rownames(R_M_miR_routine))+
  geom_smooth(method=lm)

miR223_miR_log2FC_correlation_routine <- ggplot(mapping=aes(x=sequencing_log2_FC_miRs_routine[ ,3], y=-R_M_miR_routine[ ,3]))+
  geom_point(colour="white")+
  xlab("Sequencing log2 Fold Change")+
  ylab("PCR log2 Fold Change")+
  ggtitle("miR223: miRs normalisers, routine")+
  geom_text(label=rownames(R_M_miR_routine))+
  geom_smooth(method=lm)

miR_log2FC_correlation_routine <- plot_grid(miR122_miR_log2FC_correlation_routine, miR4497_miR_log2FC_correlation_routine, miR223_miR_log2FC_correlation_routine)

#MenB correlations

miR122_miR_log2FC_correlation_MenB <- ggplot(mapping=aes(x=sequencing_log2_FC_miRs_MenB[ ,1], y=-R_M_miR_MenB[ ,1]))+
  geom_point(colour="white")+
  xlab("Sequencing log2 Fold Change")+
  ylab("PCR log2 Fold Change")+
  ggtitle("miR122: miRs normalisers, MenB")+
  geom_text(label=rownames(R_M_miR_MenB))+
  geom_smooth(method=lm)

miR4497_miR_log2FC_correlation_MenB <- ggplot(mapping=aes(x=sequencing_log2_FC_miRs_MenB[ ,2], y=-R_M_miR_MenB[ ,2]))+
  geom_point(colour="white")+
  xlab("Sequencing log2 Fold Change")+
  ylab("PCR log2 Fold Change")+
  ggtitle("miR4497: miRs normalisers, MenB")+
  geom_text(label=rownames(R_M_miR_MenB))+
  geom_smooth(method=lm)

miR223_miR_log2FC_correlation_MenB <- ggplot(mapping=aes(x=sequencing_log2_FC_miRs_MenB[ ,3], y=-R_M_miR_MenB[ ,3]))+
  geom_point(colour="white")+
  xlab("Sequencing log2 Fold Change")+
  ylab("PCR log2 Fold Change")+
  ggtitle("miR223: miRs normalisers, MenB")+
  geom_text(label=rownames(R_M_miR_MenB))+
  geom_smooth(method=lm)

miR_log2FC_correlation_MenB <- plot_grid(miR122_miR_log2FC_correlation_MenB, miR4497_miR_log2FC_correlation_MenB, miR223_miR_log2FC_correlation_MenB)


#routine/ MenB spearman rank correlations

miR122_miR_corr_routine <- cor.test(x=sequencing_log2_FC_miRs_routine[ ,1], y=-R_M_miR_routine[ ,1], method = 'spearman')
miR4497_miR_corr_routine <- cor.test(x=sequencing_log2_FC_miRs_routine[ ,2], y=-R_M_miR_routine[ ,2], method = 'spearman')
miR223_miR_corr_routine <- cor.test(x=sequencing_log2_FC_miRs_routine[ ,3], y=-R_M_miR_routine[ ,3], method = 'spearman')


miR122_miR_corr_MenB <- cor.test(x=sequencing_log2_FC_miRs_MenB[ ,1], y=-R_M_miR_MenB[ ,1], method = 'spearman')
miR4497_miR_corr_MenB <- cor.test(x=sequencing_log2_FC_miRs_MenB[ ,2], y=-R_M_miR_MenB[ ,2], method = 'spearman')
miR223_miR_corr_MenB <- cor.test(x=sequencing_log2_FC_miRs_MenB[ ,3], y=-R_M_miR_MenB[ ,3], method = 'spearman')


#TEMP CORRELATION

temp_values <- read.csv("temp correlation.csv")

temp_true <- grep("B", temp_values[ ,1])

temp_values <- temp_values[temp_true, ]

temp_order <- order(temp_values[ , 1])
temp_values <- temp_values[temp_order, ]


miR122_miR_temp_correlation <- ggplot(mapping=aes(y=temp_values[ ,2], x=-miR_DDCT[ ,1]))+
  geom_point(colour="red", shape=18, size=2)+
  ylab("Maximum Temperature, °C")+
  xlab("PCR Log2 Fold Change: miR-122")+
  geom_smooth(method=lm, colour="indianred1", fill="grey80")+
  scale_x_continuous(limit = c(-3.2, 0.8))


miR4497_miR_temp_correlation <- ggplot(mapping=aes(y=temp_values[ ,2], x=-miR_DDCT[ ,2]))+
  geom_point(colour="blue1", shape=18, size=2)+
  ylab("Maximum Temperature, °C")+
  xlab("PCR Log2 Fold Change: miR-4497")+
  geom_smooth(method=lm, colour="royalblue2", fill="grey80")

miR223_miR_temp_correlation <- ggplot(mapping=aes(y=temp_values[ ,2], x=-miR_DDCT[ ,3]))+
  geom_point(colour="green4", shape=18, size=2)+
  ylab("Maximum Temperature, °C")+
  xlab("PCR Log2 Fold Change: miR-223")+
  geom_smooth(method=lm, colour="limegreen", fill="grey80")

miR_DDCT_temp_corr <- miR_DDCT[-5, ]
temp_values_corr <- temp_values[-5, ]



miR122_miR_temp_corr <- cor.test(y=temp_values_corr[ ,2], x=-miR_DDCT_temp_corr[ ,1], method = 'spearman')
miR4497_miR_temp_corr <- cor.test(y=temp_values_corr[ ,2], x=-miR_DDCT_temp_corr[ ,2], method = 'spearman')
miR223_miR_temp_corr <- cor.test(y=temp_values_corr[ ,2], x=-miR_DDCT_temp_corr[ ,3], method = 'spearman')


#different ways of correlating fever

pre_norm_Ct_order <- order(rownames(pre_miR_normalised_Cts))
pre_miR_normalised_Cts <- pre_miR_normalised_Cts[pre_norm_Ct_order, ]

post_norm_Ct_order <- order(rownames(post_miR_normalised_Cts))
post_miR_normalised_Cts <- post_miR_normalised_Cts[post_norm_Ct_order, ]

mean_normalised_Cts <- pre_miR_normalised_Cts + post_miR_normalised_Cts
mean_normalised_Cts <- mean_normalised_Cts/2

#miR122

#pre normalised Cts vs temp
miR122_miR_pre_temp_correlation <- ggplot(mapping=aes(x=temp_values[ ,2], y=pre_miR_normalised_Cts[ ,1]))+
  geom_point(colour="black")+
  xlab("Maximum Temperature, °C")+
  ylab("PCR normalised pre Ct")+
  ggtitle("miR122: miR normalisers")+
  geom_smooth(method=lm)

miR122_miR_post_temp_correlation <- ggplot(mapping=aes(x=temp_values[ ,2], y=post_miR_normalised_Cts[ ,1]))+
  geom_point(colour="black")+
  xlab("Maximum Temperature, °C")+
  ylab("PCR normalised post Ct")+
  ggtitle("miR122: miR normalisers")+
  geom_smooth(method=lm)

miR122_miR_mean_temp_correlation <- ggplot(mapping=aes(x=temp_values[ ,2], y=mean_normalised_Cts[ ,1]))+
  geom_point(colour="black")+
  xlab("Maximum Temperature, °C")+
  ylab("PCR mean normalised Ct")+
  ggtitle("miR122: miR normalisers")+
  geom_smooth(method=lm)

#miR4497

#pre normalised Cts vs temp
miR4497_miR_pre_temp_correlation <- ggplot(mapping=aes(x=temp_values[ ,2], y=pre_miR_normalised_Cts[ ,2]))+
  geom_point(colour="black")+
  xlab("Maximum Temperature, °C")+
  ylab("PCR normalised pre Ct")+
  ggtitle("miR4497: miR normalisers")+
  geom_smooth(method=lm)

miR4497_miR_post_temp_correlation <- ggplot(mapping=aes(x=temp_values[ ,2], y=post_miR_normalised_Cts[ ,2]))+
  geom_point(colour="black")+
  xlab("Maximum Temperature, °C")+
  ylab("PCR normalised post Ct")+
  ggtitle("miR4497: miR normalisers")+
  geom_smooth(method=lm)

miR4497_miR_mean_temp_correlation <- ggplot(mapping=aes(x=temp_values[ ,2], y=mean_normalised_Cts[ ,2]))+
  geom_point(colour="black")+
  xlab("Maximum Temperature, °C")+
  ylab("PCR mean normalised Ct")+
  ggtitle("miR4497: miR normalisers")+
  geom_smooth(method=lm)

#miR223

#pre normalised Cts vs temp
miR223_miR_pre_temp_correlation <- ggplot(mapping=aes(x=temp_values[ ,2], y=pre_miR_normalised_Cts[ ,3]))+
  geom_point(colour="black")+
  xlab("Maximum Temperature, °C")+
  ylab("PCR normalised pre Ct")+
  ggtitle("miR223: miR normalisers")+
  geom_smooth(method=lm)

miR223_miR_post_temp_correlation <- ggplot(mapping=aes(x=temp_values[ ,2], y=post_miR_normalised_Cts[ ,3]))+
  geom_point(colour="black")+
  xlab("Maximum Temperature, °C")+
  ylab("PCR normalised post Ct")+
  ggtitle("miR223: miR normalisers")+
  geom_smooth(method=lm)

miR223_miR_mean_temp_correlation <- ggplot(mapping=aes(x=temp_values[ ,2], y=mean_normalised_Cts[ ,3]))+
  geom_point(colour="black")+
  xlab("Maximum Temperature, °C")+
  ylab("PCR mean normalised Ct")+
  ggtitle("miR223: miR normalisers")+
  geom_smooth(method=lm)

pre_miR_normalised_Cts_corr <- pre_miR_normalised_Cts[-5, ]
post_miR_normalised_Cts_corr <- post_miR_normalised_Cts[-5, ]
mean_normalised_Cts_corr <- mean_normalised_Cts[-5, ]

#miR122 spearman
miR122_pre_temp_corr <- cor.test(x=temp_values_corr[ ,2], y=pre_miR_normalised_Cts_corr[ ,1], method = 'spearman')
miR122_post_temp_corr <- cor.test(x=temp_values_corr[ ,2], y=post_miR_normalised_Cts_corr[ ,1], method = 'spearman')
miR122_mean_temp_corr <- cor.test(x=temp_values_corr[ ,2], y=mean_normalised_Cts_corr[ ,1], method = 'spearman')

#miR4497 spearman
miR4497_pre_temp_corr <- cor.test(x=temp_values_corr[ ,2], y=pre_miR_normalised_Cts_corr[ ,2], method = 'spearman')
miR4497_post_temp_corr <- cor.test(x=temp_values_corr[ ,2], y=post_miR_normalised_Cts_corr[ ,2], method = 'spearman')
miR4497_mean_temp_corr <- cor.test(x=temp_values_corr[ ,2], y=mean_normalised_Cts_corr[ ,2], method = 'spearman')

#miR223 spearman
miR223_pre_temp_corr <- cor.test(x=temp_values_corr[ ,2], y=pre_miR_normalised_Cts_corr[ ,3], method = 'spearman')
miR223_post_temp_corr <- cor.test(x=temp_values_corr[ ,2], y=post_miR_normalised_Cts_corr[ ,3], method = 'spearman')
miR223_mean_temp_corr <- cor.test(x=temp_values_corr[ ,2], y=mean_normalised_Cts_corr[ ,3], method = 'spearman')


