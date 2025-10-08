###############################
########installation 
###############################
devtools::install_github('Mordziarz/TAILcaller')
library(TAILcaller)
###############################
###############################
###############################
library(Rsamtools)
library(dplyr)
library(stats)
library(tidyr)
library(rlang)
library(ggplot2)
library(ComplexHeatmap)
library(ggtree)
library(circlize)
library(car)
library(dunn.test)
library(nortest)
library(rstatix)

set.seed(123) 

setwd("/dane/TAILcaller/")

################################################################################
################################ cDNA n>2 ######################################
################################################################################

###############################
#########create samples table
###############################
bamfile1 <- "CONTROL1.bam"
bamfile2 <- "CONTROL2.bam"
bamfile3 <- "CONTROL3.bam"
bamfile4 <- "CONTROL4.bam"
bamfile5 <- "CONTROL5.bam"
bamfile6 <- "HIGH1.bam"
bamfile7 <- "HIGH2.bam"
bamfile8 <- "HIGH3.bam"
bamfile9 <- "LOW1.bam"
bamfile10 <- "LOW2.bam"
bamfile11 <- "LOW3.bam"
bamfile12 <- "LOW4.bam"
bamfile13 <- "LOW5.bam"

samples_table <- data.frame(bam_path = c(bamfile1,bamfile2,bamfile3,bamfile4,bamfile5,bamfile6,bamfile7,bamfile8,bamfile9,bamfile10,
                                         bamfile11,bamfile12,bamfile13),
                            sample_name = c("CTR1","CTR2","CTR3","CTR4","CTR5","HIGH1","HIGH2","HIGH3","LOW1","LOW2","LOW3","LOW4","LOW5"),
                            group = c("CTR","CTR","CTR","CTR","CTR","HIGH","HIGH","HIGH","LOW","LOW","LOW","LOW","LOW"))
###########################
####### get info about poly(A)
############################
get_polyA_out <- get_polyA(samples_table = samples_table)
###########################
######### Linking transcripts to genes
###########################
library(rtracklayer)

gtf <- import.gff("stringtie_merge_Nano.gtf")
gtf <- as.data.frame(gtf)
get_gene_id_out <- TAILcaller::get_gene_id(polyA_table=get_polyA_out,
                                           gtf_file = gtf,
                                           transcript_column_gtf = "transcript_id",
                                           gene_column_gtf = "gene_id")
###########################
######## check duplicates
###########################
check_duplicates <- TAILcaller::polyA_duplicates(polyA_table = get_gene_id_out,
                                                 delete_duplicates = TRUE,
                                                 gene_column_gtf = "gene_id")

##########################
######### Density
##########################
density_plot <- TAILcaller::plot_density(polyA_table = get_gene_id_out,stats = "median",grouping_factor = "group")

density_plot$test
density_plot$variance
density_plot$normality
density_plot$plot

png("density_plot_3n.png", width=5, height=3, units = "in", res = 300)
density_plot$plot
dev.off()
###########################
######### Using the 'count_molecules()' function to retain genes that had more than 10 tails.
###########################

count_molecules_out <- TAILcaller::count_molecules(polyA_table = get_gene_id_out,grouping_factor="group",which_level="gene_id")
count_molecules_out <- count_molecules_out[count_molecules_out$count > 10,]

get_gene_id_out_10 <- get_gene_id_out[get_gene_id_out$gene_id %in% count_molecules_out$gene_id,]

###########################
########## PCA
###########################

count_molecules_out_pca <- TAILcaller::count_molecules(polyA_table = get_gene_id_out,grouping_factor="sample_name",which_level="gene_id")
count_molecules_out_pca <- count_molecules_out_pca[count_molecules_out_pca$count > 10,]

get_matrix_out <- TAILcaller::get_matrix(count_molecules_out = count_molecules_out_pca,
                                         grouping_factor = "sample_name",
                                         which_level = "gene_id",
                                         statistic = "median_polyA_length")

pca <- TAILcaller::PCA_polyA(get_matrix_out = get_matrix_out,samples_table = samples_table,grouping_factor = "group")

png("PCA_3n.png", width=4, height=4, units = "in", res = 300)
pca
dev.off()

###########################
######### Calculate statistics n>2
###########################

calculate_polyA_stat_n3_out <-calculate_polyA_stat_n3(polyA_table = get_gene_id_out_10,
                                                   grouping_factor = "group",
                                                   which_level = "gene_id",
                                                   padj_method="fdr")

calculate_polyA_stat_n3_out_sig <- calculate_polyA_stat_n3_out[!is.na(calculate_polyA_stat_n3_out$padj),]
calculate_polyA_stat_n3_out_sig <- calculate_polyA_stat_n3_out_sig[calculate_polyA_stat_n3_out_sig$padj < 0.05,]

write.csv2(calculate_polyA_stat_n3_out_sig,"calculate_polyA_stat_n3_out_sig.csv")


calculate_statistics_out_n3 <- TAILcaller::calculate_kruskal_polyA(polyA_table = get_gene_id_out_10,
                                                                grouping_factor = "group",
                                                                which_level = "gene_id",
                                                                padj_method="fdr")

calculate_statistics_out_n3_sig <- calculate_statistics_out_n3[!is.na(calculate_statistics_out_n3$padj),]
calculate_statistics_out_n3_sig <- calculate_statistics_out_n3_sig[calculate_statistics_out_n3_sig$padj < 0.05,]

write.csv2(calculate_statistics_out_sig,"calculate_statistics_out_n3_sig.csv")

library(ggvenn)

x <- list(`calculate_polyA_stat_n3()`=calculate_polyA_stat_n3_out_sig$gene_id,`calculate_kruskal_polyA()`=calculate_statistics_out_n3_sig$gene_id)
png(("Venn_comparison_function.png"), width=5, height=5, units = "in", res = 300)
ggvenn(x, fill_color = c("red","green","blue","yellow"), stroke_size = 0.5,set_name_size = 4,text_size = 5)
dev.off()

##########################
######### Density selected genes
##########################
polyA_table_MSTRG.961 <- get_gene_id_out_10[get_gene_id_out_10$gene_id %in% "MSTRG.961",]

density_plot_MSTRG.961 <- TAILcaller::plot_density(polyA_table = polyA_table_MSTRG.961,stats = "median",grouping_factor = "group")

density_plot_MSTRG.961$test
density_plot_MSTRG.961$variance
density_plot_MSTRG.961$normality
density_plot_MSTRG.961$plot

density_plot_MSTRG.961 <- density_plot_MSTRG.961$plot + ggtitle("MSTRG.961 ANOVA")

##############################################################################
##############################################################################
##############################################################################

polyA_table_MSTRG.20086 <- get_gene_id_out_10[get_gene_id_out_10$gene_id %in% "MSTRG.20086",]

density_plot_MSTRG.20086 <- TAILcaller::plot_density(polyA_table = polyA_table_MSTRG.20086,stats = "median",grouping_factor = "group")

density_plot_MSTRG.20086$test
density_plot_MSTRG.20086$variance
density_plot_MSTRG.20086$normality
density_plot_MSTRG.20086$plot

density_plot_MSTRG.20086 <- density_plot_MSTRG.20086$plot + ggtitle("MSTRG.20086 Kruskal-Wallis")

##############################################################################
##############################################################################
##############################################################################

polyA_table_MSTRG.19523 <- get_gene_id_out_10[get_gene_id_out_10$gene_id %in% "MSTRG.19523",]

polyA_table_MSTRG.19523 <- TAILcaller::plot_density(polyA_table = polyA_table_MSTRG.19523,stats = "median",grouping_factor = "group")

polyA_table_MSTRG.19523$test
polyA_table_MSTRG.19523$variance
polyA_table_MSTRG.19523$normality
polyA_table_MSTRG.19523$plot

polyA_table_MSTRG.19523 <- polyA_table_MSTRG.19523$plot + ggtitle("MSTRG.19523 Welch ANOVA")

############################
########## graph merge
############################
library(ggpubr)

combined_plot <- ggarrange(density_plot_MSTRG.961, density_plot_MSTRG.20086, polyA_table_MSTRG.19523,
                           common.legend = TRUE,
                           legend = "bottom",
                           ncol = 1,
                           nrow = 3)

png(("density_genes.png"), width=10, height=6, units = "in", res = 300)
combined_plot
dev.off()
###########################
######### Calculate statistics n = 2
###########################
calculate_statistics_out <-TAILcaller::calculate_statistics(polyA_table = get_gene_id_out_10,
                                                            grouping_factor = "group",
                                                            which_level = "gene_id",
                                                            control_group = "CTR",
                                                            treated_group = "HIGH")

calculate_statistics_out_sig <- calculate_statistics_out[!is.na(calculate_statistics_out$padj),]
calculate_statistics_out_sig <- calculate_statistics_out_sig[calculate_statistics_out_sig$padj < 0.05,]

write.csv2(calculate_statistics_out_sig,"calculate_statistics_out_sig.csv")

calculate_statistics_out_n2 <-TAILcaller::calculate_statistics_n2(polyA_table = get_gene_id_out_10,
                                                                  grouping_factor = "group",
                                                                  which_level = "gene_id",
                                                                  control_group = "CTR",
                                                                  treated_group = "HIGH")

calculate_statistics_out_n2_sig <- calculate_statistics_out_n2[!is.na(calculate_statistics_out_n2$padj),]
calculate_statistics_out_n2_sig <- calculate_statistics_out_n2_sig[calculate_statistics_out_n2_sig$padj < 0.05,]

write.csv2(calculate_statistics_out_n2_sig,"calculate_statistics_out_n2_sig.csv")

###########################
######### Venn n=2
###########################
x <- list(`calculate_statistics_n2()`=calculate_statistics_out_n2_sig$gene_id,`calculate_statistics()`=calculate_statistics_out_sig$gene_id)
png(("Venn_comparison_functionn=2.png"), width=5, height=5, units = "in", res = 300)
ggvenn(x, fill_color = c("blue","yellow"), stroke_size = 0.5,set_name_size = 4,text_size = 5)
dev.off()
##########################
######### Density selected genes
##########################
polyA_table_MSTRG.6021 <- get_gene_id_out_10[get_gene_id_out_10$gene_id %in% "MSTRG.6021",]
polyA_table_MSTRG.6021 <- polyA_table_MSTRG.6021[polyA_table_MSTRG.6021$group %in% c("CTR","HIGH"),]
polyA_table_MSTRG.6021 <- TAILcaller::plot_density(polyA_table = polyA_table_MSTRG.6021,stats = "median",grouping_factor = "group")

polyA_table_MSTRG.6021$test
polyA_table_MSTRG.6021$variance
polyA_table_MSTRG.6021$normality
polyA_table_MSTRG.6021$plot

polyA_table_MSTRG.6021 <- polyA_table_MSTRG.6021$plot + 
  ggtitle("MSTRG.6021 Student's t-test") + 
  scale_color_manual(values = c("CTR"="purple","HIGH"="yellow3")) 
##########################
######### Density selected genes
##########################
polyA_table_MSTRG.10967 <- get_gene_id_out_10[get_gene_id_out_10$gene_id %in% "MSTRG.10967",]
polyA_table_MSTRG.10967 <- polyA_table_MSTRG.10967[polyA_table_MSTRG.10967$group %in% c("CTR","HIGH"),]
polyA_table_MSTRG.10967 <- TAILcaller::plot_density(polyA_table = polyA_table_MSTRG.10967,stats = "median",grouping_factor = "group")

polyA_table_MSTRG.10967$test
polyA_table_MSTRG.10967$variance
polyA_table_MSTRG.10967$normality
polyA_table_MSTRG.10967$plot

polyA_table_MSTRG.10967 <- polyA_table_MSTRG.10967$plot + 
  ggtitle("MSTRG.10967 Welch's t-test") + 
  scale_color_manual(values = c("CTR"="purple","HIGH"="yellow3")) 
##########################
######### Density selected genes
##########################
polyA_table_MSTRG.3734 <- get_gene_id_out_10[get_gene_id_out_10$gene_id %in% "MSTRG.3734",]
polyA_table_MSTRG.3734 <- polyA_table_MSTRG.3734[polyA_table_MSTRG.3734$group %in% c("CTR","HIGH"),]
polyA_table_MSTRG.3734 <- TAILcaller::plot_density(polyA_table = polyA_table_MSTRG.3734,stats = "median",grouping_factor = "group")

polyA_table_MSTRG.3734$test
polyA_table_MSTRG.3734$variance
polyA_table_MSTRG.3734$normality
polyA_table_MSTRG.3734$plot

polyA_table_MSTRG.3734 <- polyA_table_MSTRG.3734$plot + 
  ggtitle("MSTRG.3734 Wilcoxon rank-sum test") + 
  scale_color_manual(values = c("CTR"="purple","HIGH"="yellow3"))
###########################
######### Volcano
###########################
volcano <- TAILcaller::volcano_polyA(calculate_statistics_out = calculate_statistics_out_n2,
                                     collapsed_color = "red",
                                     expansion_color = "green")

png("volcano.png", width=6, height=6, units = "in", res = 300)
volcano
dev.off()
###########################
######### Ma plot
###########################
ma_plot <- TAILcaller::maplot_polyA(calculate_statistics_out = calculate_statistics_out_n2,
                                    collapsed_color = "red",
                                    expansion_color = "green")

png("ma_plot.png", width=6, height=6, units = "in", res = 300)
table(ma_plot$data$PolyA_tail_length)
dev.off()
############################
########## graph merge
############################
library(ggpubr)

combined_plot_ctr_high <- ggarrange(polyA_table_MSTRG.6021, polyA_table_MSTRG.3734, polyA_table_MSTRG.10967,
                           common.legend = TRUE,
                           legend = "bottom",
                           ncol = 1,
                           nrow = 3)

png(("density_genes_ctr_high.png"), width=8, height=5, units = "in", res = 300)
combined_plot_ctr_high
dev.off()

combined_plot_ma_vol<- ggarrange(volcano, ma_plot,
                                    common.legend = TRUE,
                                    legend = "bottom",
                                    ncol = 2,
                                    nrow = 1)

png(("vol_ma.png"), width=8, height=6, units = "in", res = 300)
combined_plot_ma_vol
dev.off()

################################################################################
################################ direct RNA ####################################
################################################################################
bamfile1="L1.bam"
bamfile2="L2.bam"
bamfile3="L3.bam"
bamfile4="L4.bam"
bamfile5="W1.bam"
bamfile6="W2.bam"
bamfile7="W3.bam"
bamfile8="W4.bam"

samples_table_d2 <- data.frame(bam_path = c(bamfile1,bamfile2,bamfile3,bamfile4,bamfile5,bamfile6,bamfile7,bamfile8),
                            sample_name = c("L1","L2","L3","L4","W1","W2","W3","W4"),
                            group = c("L","L","L","L","W","W","W","W"))

###########################
####### get info about poly(A)
############################
get_polyA_out_d2 <- get_polyA(samples_table = samples_table_d2)

get_polyA_out_d2$polyA_length <- as.numeric(get_polyA_out_d2$polyA_length)
get_polyA_out_d2 <- get_polyA_out_d2[get_polyA_out_d2$polyA_length > 0,]
###########################
######### Linking transcripts to genes
###########################
library(rtracklayer)

gtf <- import.gff("stringtie_merge_Nano.gtf")
gtf <- as.data.frame(gtf)
get_gene_id_out_d2 <- TAILcaller::get_gene_id(polyA_table=get_polyA_out,
                                           gtf_file = gtf,
                                           transcript_column_gtf = "transcript_id",
                                           gene_column_gtf = "gene_id")
###########################
######## check duplicates
###########################
check_duplicates_d2 <- TAILcaller::polyA_duplicates(polyA_table = get_gene_id_out_d2,
                                                    delete_duplicates = TRUE,
                                                    gene_column_gtf = "gene_id")
check_duplicates_d2_out <- check_duplicates_d2$processed_table
##########################
######### Density
##########################
density_plot_d2 <- TAILcaller::plot_density(polyA_table = check_duplicates_d2_out,
                                            stats = "median",
                                            grouping_factor = "group")

density_plot_d2$test
density_plot_d2$variance
density_plot_d2$normality
density_plot_d2$plot

png("density_plot_2n.png", width=5, height=3, units = "in", res = 300)
density_plot_d2$plot + scale_color_manual(values = c("L"="orange","W"="blue")) 
dev.off()
###########################
######### Using the 'count_molecules()' function to retain genes that had more than 10 tails.
###########################

count_molecules_d2_out <- TAILcaller::count_molecules(polyA_table = check_duplicates_d2_out,
                                                      grouping_factor="group",
                                                      which_level="gene_id")

count_molecules_d2_out <- count_molecules_d2_out[count_molecules_d2_out$count > 10,]

get_gene_id_out_d2_10 <- check_duplicates_d2_out[check_duplicates_d2_out$gene_id %in% count_molecules_d2_out$gene_id,]

###########################
########## PCA
###########################
count_molecules_out_pca_d2 <- TAILcaller::count_molecules(polyA_table = get_gene_id_out_d2_10,
                                                          grouping_factor="sample_name",
                                                          which_level="gene_id")

get_matrix_out_d2 <- TAILcaller::get_matrix(count_molecules_out = count_molecules_out_pca_d2,
                                         grouping_factor = "sample_name",
                                         which_level = "gene_id",
                                         statistic = "median_polyA_length")

pca_d2 <- TAILcaller::PCA_polyA(get_matrix_out = get_matrix_out_d2,
                                samples_table = samples_table_d2,
                                grouping_factor = "group")

png("PCA_2n.png", width=4, height=4, units = "in", res = 300)
pca_d2 + scale_color_manual(values = c("L"="orange","W"="blue")) + geom_text(vjust=2)
dev.off()
###########################
########## Clustering/alternative visualization for density plot
###########################
polyA_heatmap_out <- TAILcaller::polyA_heatmap(polyA_table = get_gene_id_out_d2_10,
                                               grouping_factor = "sample_name",
                                               frame = 15,
                                               select="normalized",
                                               heatmap_color = "blue_green")

png("heatmap.png", width=8, height=6, units = "in", res = 300)
polyA_heatmap_out$heatmap
dev.off()
##############################
########## matrix to differentially expression analysis
##############################
count_molecules_out <- TAILcaller::count_molecules(polyA_table = get_gene_id_out_10,
                                                   grouping_factor="sample_name",
                                                   which_level="gene_id")

get_matrix_out <- TAILcaller::get_matrix(count_molecules_out = count_molecules_out,
                                         grouping_factor = "sample_name",
                                         which_level = "gene_id",
                                         statistic = "count")

write.csv(get_matrix_out,"diff_exp.csv")
