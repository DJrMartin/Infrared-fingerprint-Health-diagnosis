library(caret)
library(pROC)
library(randomForest)
load("~/Dropbox/DPSMET_2021/SET_1_2.RData")

#SET 1 AS CALIBRATION###########################################################
set.seed(12)
outliers=c(33,35,49,88)

##CLASS DEFINITION
CALIBRATION=SET_1_spectra_D2_SPLINE

y_cl = y = SET_1$Triglycerides.Hp[-outliers]
y_cl[y>30] = "Trigly_high"
y_cl[y<=30] = "Trigly_low"
y_cl=as.factor(y_cl)

CALIBRATION$Class = as.factor(y_cl)

roc_1=roc_2=roc_3=NULL
spe_1=spe_2=spe_3=NULL
sens_1=sens_2=sens_3=NULL

for (rep in (1:17)){
  #PARTITION
  inTraining <- createDataPartition(CALIBRATION$Class, p=.75, list = FALSE)
  training <- CALIBRATION[ inTraining,]
  testing  <- CALIBRATION[-inTraining,]
  #MODEL
  rf=randomForest(Class~.,data = training, mtry=200, ntree=2000)
  #PREDICTION ON VALIDATION SETS
  PRED_SET_1 = predict(rf, newdata = testing, type = "prob")
  PRED_SET_2 = predict(rf, SET_2_spectra_D2_SPLINE, type='prob')
  #ROC
  rob_1=roc(testing$Class, PRED_SET_1[,1], quiet=T)
  rob_2=roc(SET_2$Triglycerides.Hp>30,PRED_SET_2[,1], quiet=T)
  #REC
  roc_1=c(roc_1,rob_1$auc)
  roc_2=c(roc_2,rob_2$auc)
}
paste('auroc set1',median(roc_1), 
      'auroc set2',median(roc_2))
