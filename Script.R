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

set.seed(123) 

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

samples_table <- data.frame(bam_path = c(bamfile1,bamfile2,bamfile3,bamfile4,bamfile5,bamfile6,bamfile7,bamfile8),
                            sample_name = c("CTR1","CTR2","CTR3","CTR4","CTR5","HIGH1","HIGH2","HIGH3"),
                            group = c("CTR","CTR","CTR","CTR","CTR","HIGH","HIGH","HIGH"))
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
get_gene_id_out <- TAILcaller::get_gene_id(polyA_table=get_polyA_out,gtf_file = gtf,transcript_column_gtf = "transcript_id",gene_column_gtf = "gene_id")

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

png("PCA.png", width=6, height=6, units = "in", res = 300)
pca
dev.off()

###########################
######### Calculate statistics
###########################

calculate_statistics_out <-TAILcaller::calculate_statistics(polyA_table = get_gene_id_out_10,
                                                            grouping_factor = "group",
                                                            which_level = "gene_id",
                                                            control_group = "CTR",
                                                            treated_group = "HIGH")

calculate_statistics_out_sig <- calculate_statistics_out[!is.na(calculate_statistics_out$padj),]
calculate_statistics_out_sig <- calculate_statistics_out_sig[calculate_statistics_out_sig$padj < 0.05,]

###########################
######### Volcano
###########################
volcano <- TAILcaller::volcano_polyA(calculate_statistics_out = calculate_statistics_out)

png("volcano.png", width=6, height=6, units = "in", res = 300)
volcano
dev.off()

###########################
######### Ma plot
###########################
ma_plot <- TAILcaller::maplot_polyA(calculate_statistics_out = calculate_statistics_out)

png("ma_plot.png", width=6, height=6, units = "in", res = 300)
ma_plot
dev.off()

##########################
######### Density
##########################

density_plot <- TAILcaller::plot_density(polyA_table = get_polyA_out,stats = "median",grouping_column = "group")

density_plot$wilcox_test

png("density_plot.png", width=6, height=6, units = "in", res = 300)
density_plot$plot
dev.off()

###########################
########## Heatmap
###########################

res <- TAILcaller::polyA_heatmap(polyA_table = get_polyA_out,grouping_factor = "sample_name",frame = 10,select="normalized")

res$matrix

png("heatmap.png", width=6, height=6, units = "in", res = 300)
res$heatmap
dev.off()

png("tree.png", width=7.5, height=6, units = "in", res = 300)
res$tree + geom_tiplab(fontface="italic") + xlim(NA,0.9)
dev.off()
