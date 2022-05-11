# QC-backup
#The backup code to do the Quality Control
setwd('put the location information here:)')

library(dplyr)

library(data.table)

library(fastqcr)

library(ggplot2)

options(stringsAsFactors = F)

# fastQC reports
qc.dir <- list.files('QC/fastQC')

qc <- qc_aggregate(paste0('QC/fastQC/', qc.dir))

qc <- qc %>% mutate(sample = rep(gsub('-\\d+|_S.*', '', qc.dir), each=11))

qc_summ <- summary(qc)

fwrite(qc_summ, 'QC/qc_summary.csv')

qc_stat <- qc_stats(qc)

qc_stat <- as.data.frame(qc_stat) %>%
  mutate(tot.seq = as.numeric(tot.seq))

fwrite(qc_stat, 'QC/qc_stats.csv')

# Number of reads
ggsave('QC/Numreads.png', 
       ggplot(qc_stat, aes(x = sample, y = tot.seq/1e6))+
         geom_point() +
         theme_bw() +
         theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
         xlab('Samples') +
         ylab('Num of Reads (million reads)'),
       height = 4, width = 12)


# gene counts (SALMON)
genes_cnt <- readRDS('expression.rds')  # if expression matrix is ready

colnames(genes_cnt) <- gsub('-BAL-Tri-RNA.*', '', colnames(genes_cnt))

genes <- genes_cnt[rowSums(genes_cnt != 0) > 10,]

pdat <- data.frame(sample = colnames(genes)) %>%
  tidyr::separate(sample, c('id', 'visit'), remove = F)

# PCA on raw counts
pca <- prcomp(t(genes), scale. = TRUE)

ggsave('QC/PCA_rawcounts.png',
       ggbiplot::ggbiplot(pca, obs.scale = 1, var.scale = 1,
                          ellipse = TRUE, varname.size=TRUE, pc.biplot=FALSE, var.axes=FALSE,
                          circle = FALSE,labels = colnames(genes),
                          group = pdat$visit) + theme_bw()+
         theme(legend.position = 'bottom')+
         scale_color_manual(values = c("#e41a1c", "#377eb8"))+
         theme(axis.title=element_text(face="bold",size="14"),
               axis.text=element_text(face="bold",size="14"),
               legend.text=element_text(face="bold",size="11"),
               legend.title=element_blank(),
               panel.grid.major = element_blank(),
               panel.grid.minor = element_blank())+
         guides(color=guide_legend(nrow=1,byrow=TRUE)),
       height = 8, width = 8)

# boxplot on raw counts
png('QC/log10_rawcounts.png', width = 1000, height = 800)

boxplot(log10(genes+1), las=2)

dev.off()

# library size
plotdat = data.frame(Sum=colSums(genes_cnt), sampl=colnames(genes_cnt))

ggsave('QC/Libsize.png', 
       ggplot(plotdat, aes(x=sampl, y=Sum/1e6))+
         geom_point()+
         theme_bw()+
         theme(axis.text.x = element_text(angle = 90, hjust = 1))+
         xlab('Samples')+
         ylab('Lib.Size (million reads)'),
       height=4, width=12)

# PCA against seq date
sdate <- readLines('info/DISARM_BAL_seq_date.txt')

sdate <- data.frame(sample = sdate[seq(1, 239, by = 2)], 
                    date = sdate[seq(2, 240, by = 2)])

sdate <- filter(sdate, grepl('DIS', sample)) %>%
  mutate(sample = gsub('-BAL-Tri-RNA.*', '', sample))

genes <- genes[, sdate$sample]

pca <- prcomp(t(genes), scale. = TRUE)

ggsave('QC/PCA_rawcounts_seq_date.png',
       ggbiplot::ggbiplot(pca, obs.scale = 1, var.scale = 1,
                          ellipse = TRUE, varname.size=TRUE, pc.biplot=FALSE, var.axes=FALSE,
                          circle = FALSE,labels = colnames(genes),
                          group = sdate$date) + theme_bw()+
         theme(legend.position = 'bottom')+
         theme(axis.title=element_text(face="bold",size="14"),
               axis.text=element_text(face="bold",size="14"),
               legend.text=element_text(face="bold",size="11"),
               legend.title=element_blank(),
               panel.grid.major = element_blank(),
               panel.grid.minor = element_blank())+
         guides(color=guide_legend(nrow=1,byrow=TRUE)),
       height = 8, width = 8)

