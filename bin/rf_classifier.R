# RANDOM FOREST CLASSIFIER

# load packages
library(randomForest);

# read arguments
write("# Reading arguments\n",stderr());
arguments <- commandArgs(TRUE);
RF     <- arguments[1];
INPUT  <- arguments[2];
OUTPUT <- arguments[3];

# load random forest
write("# Loading Random Forest\n",stderr());
load(RF);

# read table
write("# Reading data\n",stderr());
data <- read.table(file=INPUT, header=T, sep="\t", comment.char = "");

# Remove self-homologous-interactions
data <- data[as.character(data$HOM_1) != as.character(data$HOM_2),]

# Discretization!
write("# Discretizing data\n",stderr());
# REMOVE NA's before discretization!
data <- data[!is.na(data$PATH_LENGTH),]

# Can't use MDLP -> too slow
# Drop unnecessary variables and Convert to factors some of the variables
data.clean <- data[ , !(names(data) %in% c("TRANS_1", "TRANS_2", "HOM_1", "HOM_2"))]
data.clean$BLAST_BRH_1 <- factor(data.clean$BLAST_BRH_1, labels=c("no", "yes"))
data.clean$BLAST_BRH_2 <- factor(data.clean$BLAST_BRH_2, labels=c("no", "yes"))
data.clean$NOG_BRH_1   <- factor(data.clean$NOG_BRH_1, labels=c("no", "yes"))
data.clean$NOG_BRH_2   <- factor(data.clean$NOG_BRH_2, labels=c("no", "yes"))
data.clean$PFAM_BRH_1  <- factor(data.clean$PFAM_BRH_1, labels=c("no", "yes"))
data.clean$PFAM_BRH_2  <- factor(data.clean$PFAM_BRH_2, labels=c("no", "yes"))
data.clean$PATH_LENGTH <- factor(data.clean$PATH_LENGTH)

# Evalue logarithm
data.clean$BLAST_EVAL_1 <- log10(data.clean$BLAST_EVAL_1)
data.clean$BLAST_EVAL_2 <- log10(data.clean$BLAST_EVAL_2)
data.clean$NOG_EVAL_1   <- log10(data.clean$NOG_EVAL_1)
data.clean$NOG_EVAL_2   <- log10(data.clean$NOG_EVAL_2)

write(paste("Before normalization:\n", summary(data.clean)),stderr());
# Discretize evalues manually!
# 7 LEVELS!
data.clean$BLAST_EVAL_1 <- cut(data.clean$BLAST_EVAL_1, c(-Inf, -200, -150, -100, -50, -20, 0, 1), include.lowest=T)
data.clean$BLAST_EVAL_2 <- cut(data.clean$BLAST_EVAL_2, c(-Inf, -200, -150, -100, -50, -20, 0, 1), include.lowest=T)
data.clean$NOG_EVAL_1   <- cut(data.clean$NOG_EVAL_1,   c(-Inf, -200, -150, -100, -50, -20, 0, 1), include.lowest=T)
data.clean$NOG_EVAL_2   <- cut(data.clean$NOG_EVAL_2,   c(-Inf, -200, -150, -100, -50, -20, 0, 1), include.lowest=T)

# Discretize coverage manually
data.clean$BLAST_COV_1 <- cut(data.clean$BLAST_COV_1, c(0, 25, 50, 75, Inf), include.lowest=T)
data.clean$BLAST_COV_2 <- cut(data.clean$BLAST_COV_2, c(0, 25, 50, 75, Inf), include.lowest=T)
data.clean$PFAM_COV_1  <- cut(data.clean$PFAM_COV_1,  c(0, 25, 50, 75, Inf), include.lowest=T)
data.clean$PFAM_COV_2  <- cut(data.clean$PFAM_COV_2,  c(0, 25, 50, 75, Inf), include.lowest=T)

# Discretize DOMAIN INT SCORE
data.clean$DOMAIN_INT_SCORE <- cut(data.clean$DOMAIN_INT_SCORE, c(0, 1, 2, 3, 4, 5, Inf), right=F)

# Discretize NTO
data.clean$MOLFUN_NTO  <- cut(data.clean$MOLFUN_NTO,  c(-1, 0, 0.25, 0.5, 0.75, 1, Inf), include.lowest=T, right=F)
data.clean$BIOPROC_NTO <- cut(data.clean$BIOPROC_NTO, c(-1, 0, 0.25, 0.5, 0.75, 1, Inf), include.lowest=T, right=F)
data.clean$CELLCOM_NTO <- cut(data.clean$CELLCOM_NTO, c(-1, 0, 0.25, 0.5, 0.75, 1, Inf), include.lowest=T, right=F)

levels(data.clean$PATH_LENGTH) <- c("1", "2", "-1",  "3", "4", "5", "6")
write(paste("After normalization:\n", summary(data.clean)),stderr());

# PREDICTION
write(paste("# Predicting ", nrow(data.clean), " interactions...", sep=""), stderr());
p <- predict(train.rf4, data.clean, type="prob");
newdata <- cbind(data, p);
write("done\n",stderr());

write("# Writing output\n",stderr());
write.table(newdata, file=OUTPUT, sep="\t", quote=F, row.names=F);
