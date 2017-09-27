#!/usr/bin/R
library(ggplot2);
arguments <- commandArgs(TRUE);
INPUT  <- arguments[1];
OUTPUT <- arguments[2];
pfam_score <- read.table(file=INPUT, header=TRUE, sep="\t",comment.char = "");
gplot <- ggplot(pfam_score) +
            geom_bar(
                aes(x=SCORE),
                binwidth=.05
            ) +
            theme_bw() +
            scale_x_log10() +
            ylab("PFAM BRHs\n") +
            xlab("Meta-alignment Score");


ggsave(file=OUTPUT, gplot);
