installed = installed.packages()
if (!("caret" %in% installed)) { install.packages("caret") }
if (!("pROC" %in% installed)) { install.packages("pROC") }
if (!("randomForest" %in% installed)) { install.packages("randomForest") }

library(caret)
library(pROC)
library(randomForest)

# DATA IMPORTATION
load("~/Dropbox/DPSMET_2021/SET_1_2.RData")

# SET 1 USED AS CALIBRATION SET ###########################################################
set.seed(12)
outliers=c(33,35,49,88)

## CLASS DEFINITION < or > of 30mg of triglycerides per g of liver
CALIBRATION = SET_1_spectra_D2_SPLINE

y_cl = y = MIR.1$Triglycerides.Hp[-outliers]
y_cl[y>30] = "Trigly_high"
y_cl[y<=30] = "Trigly_low"
y_cl=as.factor(y_cl)

CALIBRATION$Class = as.factor(y_cl)

roc_1=roc_2=NULL

# Classification using RANDOM FOREST ALGORITHM
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
  rob_1 = roc(testing$Class, PRED_SET_1[,1], quiet=T)
  rob_2 = roc(MIR.2$Triglycerides.Hp>30,PRED_SET_2[,1], quiet=T)
  #REC
  roc_1 = c(roc_1,rob_1$auc)
  roc_2 = c(roc_2,rob_2$auc)
}
paste('auroc set1',median(roc_1), 
      'auroc set2',median(roc_2))

# PLOT OF THE MOST 15 IMPORTANT VARIABLES
row_names_all = NULL
for (i in c(7,5,2)){
  row_names = freq_set1[fingerprint_set1[seq(1,length(fingerprint_set1),i)]]
  row_names_all = c(row_names_all,row_names)
}
colnames(SET_1_spectra_D2_SPLINE) = row_names_all[1:dim(SET_1_spectra_D2_SPLINE)[2]]
colnames(SET_2_spectra_D2_SPLINE) = row_names_all[1:dim(SET_2_spectra_D2_SPLINE)[2]]
rownames(rf$importance) = row_names_all[1:length(rf$importance)]

require(corrplot)
# CORRELATION BETWEEN ABSORBANCES AT SELECTED WAVENUMBERS AND OUTPUT
fiveteen_variables=order(rf$importance, decreasing = TRUE)[1:10]
Imp_var_data = cbind(SET_1_spectra_D2_SPLINE[,fiveteen_variables],MIR.1$Triglycerides.Hp[-outliers],MIR.1$Steatosis[-outliers])
colnames(Imp_var_data) = c(colnames(SET_1_spectra_D2_SPLINE[,fiveteen_variables]),'Hep. Trigly.','Steatosis')
Imp_var_data_val=cbind(SET_2_spectra_D2_SPLINE[,fiveteen_variables],MIR.2$Triglycerides.Hp ,MIR.2$Steatosis)
colnames(Imp_var_data_val)=c(colnames(SET_2_spectra_D2_SPLINE[,fiveteen_variables]),'Hep. Trigly.','Steatosis')
r=rbind(Imp_var_data[, ], Imp_var_data_val[,])
C.Imp = data.frame(cor(rbind(Imp_var_data[, ], Imp_var_data_val[,])))
new=as.matrix(C.Imp$Hep..Trigly.[c(1:10)])
rownames(new)=c(colnames(Imp_var_data)[1:10])
colnames(new)='HTG'

## CORRELOGRAMME PLOT
par(fig=c(0,0.6,0.02,0.98))
library(plot.matrix)
plot(as.matrix(new), main='', xlab='', ylab='', yaxt='n',xaxt = "n" ,cex.main=0.7,
     border=NA, text.cell=list(cex=0.7), digits=2, 
     col=brewer.pal(n = 9, name = 'BuGn'), key=NULL)
## PLOT OF R SQUARED OF EACH WAVENUMBERS AND TRIGLY
R_squared= R_Squared2=p_value=NULL
for(i in 1:10){
  R=summary(lm(r$`Hep. Trigly.`~r[,i]))
  R_squared=c(R_squared,R$r.squared)
  p_value=c(p_value,R$coefficients[2,4])
}
par(fig=c(0.40,1,0,1), new=TRUE)
barplot(rev(R_squared),col="grey", horiz = TRUE, main='', cex.main=0.7, density=20)
