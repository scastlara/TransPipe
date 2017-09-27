#!/usr/bin/R
library(ggplot2);
library(scales);

orflength <- read.table(file="INFOSEQ/longestorf.infoseq", header=T, comment.char = "");
trlength <- read.table(file="INFOSEQ/transcriptome.infoseq", header=T, comment.char = "");

merged <- merge(trlength, orflength, by="Name");
colnames(merged) <- c("NAME", "TRANSCRIPT", "ORF");
merged$TRANSCRIPT <- merged$TRANSCRIPT / 3;

# BREAKS
all_lengths <- c(merged$TRANSCRIPT, merged$TORF);
breaks <- pretty(log10(all_lengths), n=5);
breaks2 <- as.integer(10^breaks);
RoundUp <- function(from,to) ceiling(from/to)*to
breaks2 <- RoundUp(breaks2, 10);


plot <- ggplot(merged, aes(x=TRANSCRIPT, y=ORF)) +
    geom_point(alpha = 0.5, color="#EBD170") +
    xlab("\nTranscript Length (number of codons)") +
    ylab("Protein Length (number of aa)\n") +
    stat_density2d(color="#002640") +
    scale_x_log10(breaks = breaks2) + scale_y_log10(breaks = breaks2) +
    stat_function(fun=function(x)x+log10(1), geom="line", colour="#4d4d4d", alpha=0.7, linetype="dashed") +
    stat_function(fun=function(x)x+log10(0.75), geom="line", colour="#4d4d4d",alpha=0.7, linetype="dashed") +
    stat_function(fun=function(x)x+log10(0.5), geom="line", colour="#4d4d4d",alpha=0.7, linetype="dashed") +
    stat_function(fun=function(x)x+log10(0.25), geom="line", colour="#4d4d4d",alpha=0.7, linetype="dashed") +
    stat_function(fun=function(x)x+log10(0.1), geom="line", colour="#4d4d4d",alpha=0.7, linetype="dashed") +
    annotate("text", label=c("100%","75%","50%","25%", "10%"), x=25000, y=c(25000,18750,12500,6250, 2500),
             angle=45, hjust=c(0.9,1.2,1.2,1.2,1.2), vjust=c(-0.1,1.1,1.1,1.1,1.1), size=3.5) +
    theme_bw()

ggsave(file="PLOTS/LONGEST_ORFS.svg", plot);
