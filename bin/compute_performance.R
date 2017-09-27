library(ROCR);
compute.AUC <- function(rf,class) {

    predictions=as.vector(rf$votes[,2])
    pred=prediction(predictions,class)

    perf_AUC=performance(pred,"auc") #Calculate the AUC value
    AUC=perf_AUC@y.values[[1]]

    perf_ROC=performance(pred,"tpr","fpr") #plot the actual ROC curve

    plot(perf_ROC, main="ROC plot")
    text(0.5,0.5,paste("AUC = ",format(AUC, digits=5, scientific=FALSE)))
    abline(0,1,lty=4)

    return(AUC)
}

compute.performance <- function(rf, class) {
    TP <- rf$confusion[,2][2]
    TN <- rf$confusion[,1][1]
    FP <- rf$confusion[,2][1]
    FN <- rf$confusion[,1][2]

    precision   <- TP / (TP + FP)
    recall      <- TP / (TP + FN)
    specificity <- TN / (TN + FP)
    f.measure   <- 2*(precision * recall) / (precision + recall)
    accuracy    <- (TP + TN) / (TP + TN + FP + FN)

    results <- data.frame(precision=precision, recall=recall, specificity=specificity, f.measure=f.measure, accuracy=accuracy)
    results$AUC <- compute.AUC(rf, class)
    return(results)
}

