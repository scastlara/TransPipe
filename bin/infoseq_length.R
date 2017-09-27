#!/usr/bin/R

library(ggplot2);

# READ DATA
subj_df  <- read.table(file="INFOSEQ/subject.infoseq", header=TRUE);
trans_df <- read.table(file="INFOSEQ/transcriptome.infoseq", header=TRUE, comment.char = "");

# PUT SOURCE ON DATA FRAMES
trans_df$SOURCE <- "TRANSCRIPTOME";
subj_df$SOURCE  <- "SUBJECT";

# CORRECT TRANSCRIPT LENGTH
trans_df$Length <- trans_df$Length / 3;

# JOIN DATA FRAMES
DF <- rbind(subj_df,trans_df)
DF$SOURCE <- as.factor(DF$SOURCE);

# BREAKS
quantiles      <- quantile(DF$Length);
vquantiles     <- as.numeric(quantiles);
vquantiles     <- c(vquantiles[1], vquantiles[3], vquantiles[5]);
logquantiles   <- log10(vquantiles);
interquantiles <- c((logquantiles[1] + logquantiles[2])/2, (logquantiles[2] + logquantiles[3])/2);
interquantiles <- c(10^interquantiles[1], 10^interquantiles[2]);

vbreaks <- c(vquantiles[1], interquantiles[1], vquantiles[2], interquantiles[2], vquantiles[3]);
vbreaks <- as.integer(vbreaks);

# PLOT THE THING
ggplot(DF, aes(x=SOURCE, y=Length, fill=SOURCE)) +
    geom_violin(alpha=0.5)   +
    geom_boxplot(width=0.15)  +
    theme_bw()      +
    scale_y_log10(breaks=c(vbreaks)) +
    theme(legend.position="none") +
    xlab("") + ylab("Length (aa or codons)\n")

ggsave(file="PLOTS/SEQUENCE_LENGTH.svg");
