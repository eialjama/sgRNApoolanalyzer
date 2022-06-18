library(tidyverse)
library(readxl)
library(ggplot2)
library(unikn)
library(readxl)
library(plyr)
library(dplyr)

setwd("C:/Users/matthew.woods/Desktop/Ion Torrent Read Counter and Aligner/Density Plot")

EaC <- read_excel("Exact alignment Counts.xlsx")
pat <- "(.*)-(.*)"
colnames(EaC) <- c('matching_Sequence', 'Count', 'Crispr_ID', 'Sequence')


EAG <- transform(EaC, Crispr_ID = sub(pat, "\\1", Crispr_ID))
EAGC <- ddply(EAG,"Crispr_ID", numcolwise(sum))

EAG <- EAG %>%
  arrange(Count)

EAGC <- EAGC %>%
  arrange(Count)


################################################################################


log_graph_gRNA <- log(EAG$Count, base = exp(2))             #by gRNA
log_graph_gRNA <- na.omit(log_graph_gRNA)                    #by gRNA


barplot(log_graph_gRNA,
        ylim = c(0,6),
        ylab = 'Base Read Depth',
        xlab =  nrow(EaC))



logDensity_gRNA <- density(log_graph_gRNA)                   #by gRNA


LogquantHigh95_gRNA <- quantile(log_graph_gRNA, probs = .95)  # by gRNA
LogquantLow05_gRNA <- quantile(log_graph_gRNA, probs = .05)


fold_95_gRNA <- LogquantHigh95_gRNA/LogquantLow05_gRNA  # by gRNA





################################################################################
#log_graph_gRNA_IO <- log(CrisprID$Ion, base = exp(2))       #by gRNA
#log_graph_gRNA_IO <- na.omit(log_graph_gRNA_IO)                    #by gRNA

#logDensity_gRNA_IO <- density(log_graph_gRNA_IO)                   #by gRNA


#LogquantHigh95_gRNA_IO <- quantile(log_graph_gRNA_IO, probs = .95)  # by gRNA
#LogquantLow05_gRNA_IO <- quantile(log_graph_gRNA_IO, probs = .05)


#fold_95_gRNA_IO <- LogquantHigh95_gRNA_IO/LogquantLow05_gRNA_IO  # by gRNA


#

# by gRNA


plot(logDensity_gRNA, 
     main= 'Prevalancy of Total Reads Scanned for Entire Library by Crispr-ID', 
     lwd = 3, 
     col = "blue", 
     xlim = c(0,10), 
     ylim = c(0, 3), 
     yaxs = "i", 
     xaxs = "i", 
     ylab = expression(paste("Frequency = (Probability of any Total Read) ", Log[2])), 
     xlab = expression(paste(" Normalized Total Reads", " = (Total Reads)", Log[2],)), 
     cex.axis = 1.1, 
     cex.lab = 1.05,
     cex.main= 1.2)


#lines(logDensity_gRNA, col="red",lwd=3)

  
annotate("text", x = 5, y = 5, parse = TRUE,
         label = "x=Rlog[2]")

par(mar = c(5,5,4,5))

abline(v = LogquantHigh95_gRNA, col = "black") # by gRNA
abline(v = LogquantLow05_gRNA, col = "black")  # by gRNA

#abline(v = LogquantHigh95_gRNA_IL, col = "purple") # by gRNA
#abline(v = LogquantLow05_gRNA_IL, col = "purple")  # by gRNA





################################################################################

log_graph_gRNA <- log(EAGC$Count, base = exp(2))  #by gene
log_graph_gRNA <- na.omit(log_graph_gRNA)                    #by gene

logDensity_gene <- density(log_graph_gRNA)                   #by gene

LogquantHigh99_gene <- quantile(log_graph_gRNA, probs = .99)  # by gene
LogquantLow01_gene <- quantile(log_graph_gRNA, probs = .1)

fold_99_gene <- LogquantHigh99_gene/LogquantLow01_gene  # by gene



################################################################################
#log_graph_gRNA_IL <- log(Gene$Illumina, base = exp(2))  #by gene
#log_graph_gRNA_IL <- na.omit(log_graph_gRNA_IL)                    #by gene

#logDensity_gene_IL <- density(log_graph_gRNA_IL)                   #by gene

#LogquantHigh99_gene_IL <- quantile(log_graph_gRNA_IL, probs = .99)  # by gene
#LogquantLow01_gene_IL <- quantile(log_graph_gRNA_IL, probs = .1)

#fold_99_gene_IL <- LogquantHigh99_gene_IL/LogquantLow01_gene_IL  # by gene


################################################################################
# by gene
plot(logDensity_gene, 
     main= 'Prevalancy of Total Reads Scanned for Entire Library by Gene', 
     lwd = 3, 
     col = "blue", 
     xlim = c(0,10), 
     ylim = c(0, 5), 
     yaxs = "i", 
     xaxs = "i", 
     ylab = expression(paste("Frequency = (Probability of any Total Read) ", Log[2])), 
     xlab = expression(paste(" Normalized Total Reads", " = (Total Reads)", Log[2],)), 
     cex.axis = 1.1, 
     cex.lab = 1.05,
     cex.main= 1.2)

#lines(logDensity_gene_IL, col="red",lwd=3)

annotate("text", x = 5, y = 5, parse = TRUE,
         label = "x=Rlog[2]")
par(mar = c(5,5,4,5))

abline(v = LogquantHigh99_gene, col = "black") # by gRNA
abline(v = LogquantLow01_gene, col = "black")  # by gRNA

#abline(v = LogquantHigh99_gene_IL, col = "purple") # by gRNA
#abline(v = LogquantLow01_gene_IL, col = "purple")  # by gRNA






