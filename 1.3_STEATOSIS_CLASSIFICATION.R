library(caret)
library(pROC)
library(randomForest)

#DATA IMPORTATION
load("~/Dropbox/DPSMET_2021/SET_1_2.RData")

#SET 1 AS CALIBRATION###########################################################
set.seed(1)
outliers=c(33,35,49,88)
na=which(is.na(MIR.1$Steatosis))

##CLASS DEFINITION
CALIBRATION=SET_1_spectra_D2_SPLINE[-na,]
y_cl = y = as.numeric((MIR.1$Steatosis[-c(na,outliers)]))*100
y_cl[y>10] = "Steatosis_high"
y_cl[y<=10] = "Steatosis_low"
y_cl=as.factor(y_cl)

CALIBRATION$Class = as.factor(y_cl)

roc_1=roc_2=NULL

for (rep in (1:20)){
  #PARTITION
  inTraining <- createDataPartition(CALIBRATION$Class, p=0.75, list = FALSE)
  training <- CALIBRATION[ inTraining,]
  testing  <- CALIBRATION[-inTraining,]
  #MODEL
  rf=randomForest(Class~.,data = training, mtry=200, ntree=2000)
  #PREDICTION ON VALIDATION SETS
  PRED_SET_1 = predict(rf, newdata = testing, type = "prob")
  PRED_SET_2 = predict(rf, SET_2_spectra_D2_SPLINE, type='prob')
  #ROC
  rob_1=roc(testing$Class, PRED_SET_1[,1], quiet=T)
  rob_2=roc(MIR.2$Steatosis>10,PRED_SET_2[,1], quiet=T)
  #ROC
  roc_1=c(roc_1,rob_1$auc)
  roc_2=c(roc_2,rob_2$auc)
}
paste('auroc set1',mean(roc_1), 
      'auroc set2',mean(roc_2))
