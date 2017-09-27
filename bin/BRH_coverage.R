library(ggplot2);
library(svglite);
library(ggpubr);
arguments <- commandArgs(TRUE);
INPUT   <- arguments[1];
OUTPUT  <- arguments[2];
PROGRAM <- arguments[3];
QUERY   <- arguments[4];
SUBJECT <- arguments[5];


BRHcov <- read.table(file=INPUT, header=F, comment.char = "");

BRHcov <- read.table(file=INPUT, header=TRUE, comment.char = "");

scatt<-ggplot(BRHcov, aes(x=TRANS_COV, y=SUBJ_COV) ) +
           geom_point(color="#ef8f2c", alpha=0.5) +
           xlab(paste("\n", QUERY, " Coverage\n")) +
           ylab(paste(SUBJECT, " Coverage\n")) +
           xlim(0,100) + ylim(0,100) +
           theme_bw()  +
           theme(legend.position = "none") +  stat_density2d(color="black", alpha=0.8) +
           geom_abline(aes(slope=1, intercept=0), linetype="dashed");


trans<-ggplot(BRHcov, aes(x=TRANS_COV)) +
        geom_bar(fill="#3d9fac", alpha=0.8, binwidth=1, position="identity") +
        xlab(paste("\n", QUERY, " Coverage\n")) + ylab("Count\n") + xlim(0,100) +
        geom_vline(xintercept= median(BRHcov$TRANS_COV), linetype="dashed") +
        theme_bw();

subj<-ggplot(BRHcov, aes(x=SUBJ_COV)) +
        geom_bar(fill="#ef5030", alpha=0.8, binwidth=1, position="identity") +
        xlab(paste(SUBJECT, " Coverage\n")) + ylab("\nCount") + xlim(0,100) +
        geom_vline(xintercept= median(BRHcov$SUBJ_COV), linetype="dashed") +
        theme_bw() + coord_flip();

title <- ggparagraph(text=paste("BRHs ", PROGRAM, "\nn = ", nrow(BRHcov)),
                                       color = "black", face = "bold", size = 12, );
arranged <- ggarrange(trans + rremove("x.title"), title, scatt, subj + rremove("y.title"), ncol=2, nrow=2, align="hv");
ggexport(arranged, filename=OUTPUT, width=500, height=500);
    #png(filename=OUTPUT);
    #arrangeGrob(trans + rremove("x.title"), junk, scatt, subj + rremove("y.title"), ncol=2);
    #dev.off();
    #gplots <- arrangeGrob(trans, junk, scatt, subj, ncol=2);
    #class(gplots) <- c("arrange","ggplot", class(gplots));
    #print.arrange <- function(x) grid.draw(x);
    #ggsave(file=OUTPUT, gplots);
#ggsave(file=OUTPUT, gplots);
