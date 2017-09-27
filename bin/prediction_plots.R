library(ggplot2)

arguments <- commandArgs(TRUE);
FILE <- arguments[1];
data <- read.table(FILE, header=T, sep="\t", comment.char = "");


# Drop unnecessary variables and Convert to factors some of the variables
data <- data[ , !(names(data) %in% c("TRANS_1", "TRANS_2", "HOM_1", "HOM_2"))]
data$BLAST_BRH_1 <- factor(data$BLAST_BRH_1, labels=c("no", "yes"))
data$BLAST_BRH_2 <- factor(data$BLAST_BRH_2, labels=c("no", "yes"))
data$NOG_BRH_1   <- factor(data$NOG_BRH_1, labels=c("no", "yes"))
data$NOG_BRH_2   <- factor(data$NOG_BRH_2, labels=c("no", "yes"))
data$PFAM_BRH_1  <- factor(data$PFAM_BRH_1, labels=c("no", "yes"))
data$PFAM_BRH_2  <- factor(data$PFAM_BRH_2, labels=c("no", "yes"))
data$PATH_LENGTH <- factor(data$PATH_LENGTH)

# Evalue logarithm
data$BLAST_EVAL_1 <- log10(data$BLAST_EVAL_1)
data$BLAST_EVAL_2 <- log10(data$BLAST_EVAL_2)
data$NOG_EVAL_1   <- log10(data$NOG_EVAL_1)
data$NOG_EVAL_2   <- log10(data$NOG_EVAL_2)

write(paste("Before normalization:\n", summary(data)),stderr());
# Discretize evalues manually!
# 7 LEVELS!
data$BLAST_EVAL_1 <- cut(data$BLAST_EVAL_1, c(-Inf, -200, -150, -100, -50, -20, 0, 1), include.lowest=T)
data$BLAST_EVAL_2 <- cut(data$BLAST_EVAL_2, c(-Inf, -200, -150, -100, -50, -20, 0, 1), include.lowest=T)
data$NOG_EVAL_1   <- cut(data$NOG_EVAL_1,   c(-Inf, -200, -150, -100, -50, -20, 0, 1), include.lowest=T)
data$NOG_EVAL_2   <- cut(data$NOG_EVAL_2,   c(-Inf, -200, -150, -100, -50, -20, 0, 1), include.lowest=T)

# Discretize coverage manually
data$BLAST_COV_1 <- cut(data$BLAST_COV_1, c(0, 25, 50, 75, Inf), include.lowest=T)
data$BLAST_COV_2 <- cut(data$BLAST_COV_2, c(0, 25, 50, 75, Inf), include.lowest=T)
data$PFAM_COV_1  <- cut(data$PFAM_COV_1,  c(0, 25, 50, 75, Inf), include.lowest=T)
data$PFAM_COV_2  <- cut(data$PFAM_COV_2,  c(0, 25, 50, 75, Inf), include.lowest=T)

# Discretize DOMAIN INT SCORE
data$DOMAIN_INT_SCORE <- cut(data$DOMAIN_INT_SCORE, c(0, 1, 2, 3, 4, 5, Inf), right=F)

# Discretize NTO
data$MOLFUN_NTO  <- cut(data$MOLFUN_NTO,  c(-1, 0, 0.25, 0.5, 0.75, 1, Inf), include.lowest=T, right=F)
data$BIOPROC_NTO <- cut(data$BIOPROC_NTO, c(-1, 0, 0.25, 0.5, 0.75, 1, Inf), include.lowest=T, right=F)
data$CELLCOM_NTO <- cut(data$CELLCOM_NTO, c(-1, 0, 0.25, 0.5, 0.75, 1, Inf), include.lowest=T, right=F)

pred <- ifelse(data$TRUE. >= 0.6, "TRUE", "FALSE")
data$pred <- pred

for (i in c("MOLFUN_NTO", "BIOPROC_NTO", "CELLCOM_NTO", "MOLFUN_NTO", "DOMAIN_INT_SCORE", "PATH_LENGTH")) {
    gplot <- ggplot(data) +
    geom_bar(aes_string(x="pred", fill=i), position="fill") +
    theme_bw() +
    ylab("Proportion of labelled Pairs\n") +
    xlab("") +
    guides(fill=guide_legend(title=i)) +
    scale_x_discrete(labels=c("FALSE", "TRUE"))

    ggsave(gplot, file=paste("PLOTS/PREDICTIONS_", i, ".png", sep=""))
}
