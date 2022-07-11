installed = installed.packages()
if (!("caret" %in% installed)) { install.packages("caret") }
if (!("pROC" %in% installed)) { install.packages("pROC") }
if (!("randomForest" %in% installed)) { install.packages("randomForest") }
if (!("corrplot" %in% installed)) { install.packages("corrplot") }

library(caret)
library(pROC)
library(randomForest)

#LOAD META DATA AND SPLINES COMPRESSION
setwd(dir="~/Desktop/DATA_HEALTH_DIAGNOSIS")
dir()
##METADATA
MIR.1 = read.csv("Blood_smear_bruts_SET1.csv", sep=';' )
MIR.2 = read.csv("Blood_smear_bruts_SET2.csv", sep=';')
##B_SPLINES COMPRESSION
load('SET_1_2_BLOOD_SMEAR.RData')

#DATA OF INFLAMMATION MARKERS IN LIVER
INFLAM_SET1 = data.frame(scale(MIR.1[, c(19,20,22)]))
INFLAM_SET2 = data.frame(scale(MIR.2[, c(14:16)]))

na = which(is.na(INFLAM_SET1$NRF2))
SET_1_spectra_D2_SPLINE = SET_1_spectra_D2_SPLINE[-na,]
y = as.numeric(as.character(MIR.1$group_inflam))[-na]

set.seed(12) # seed is fixed for repetability
y_cl = y
y_cl[y==1] = "inflammation_high"
y_cl[y==2] = "inflammation_low"
y_cl=as.factor(y_cl)
SET_1_spectra_D2_SPLINE$Class = y_cl

pred_all = NULL
roc_all = roc_all2 = NULL

#DEFINITION OF INFLAMMATION CLASS FOR SET 2
group_set2 = (INFLAM_SET2$A2M.1>-.18)&(INFLAM_SET2$CRP.1>-.3)
y_cl = group_set2
group_set2[y_cl==TRUE] = "inflammation_high"
group_set2[y_cl==FALSE]  = "inflammation_low"
group_set2 = as.factor(group_set2)

for (rep in (1:10)){
  inTraining <- createDataPartition(SET_1_spectra_D2_SPLINE$Class, p=.75, list = FALSE)
  training <- SET_1_spectra_D2_SPLINE[ inTraining,]
  testing  <- SET_1_spectra_D2_SPLINE[-inTraining,]
  rf=randomForest(Class~.,data = training, mtry=30, ntree=2000)
  # Set 1 (validation part)
  pred = predict(rf, newdata = testing, type = "prob")
  rob=roc(testing$Class, pred[,1],quiet=TRUE)
  roc_all=c(roc_all,rob$auc)
  # Set 2
  pred_val = predict(rf, newdata = SET_2_spectra_D2_SPLINE, type='prob')
  rob2=roc(group_set2, pred_val[,1], quiet=T)
  roc_all2=c(roc_all2,rob2$auc)
}
mean(roc_all)
mean(roc_all2)

##PLOT OF THE MOST 15 IMPORTANT VARIABLES
row_names_all = NULL
for (i in spline_succession){
  row_names = freq_set1[fingerprint_set1[seq(1,length(fingerprint_set1),i)]]
  row_names_all = c(row_names_all,row_names)
}
colnames(SET_1_spectra_D2_SPLINE) = row_names_all[1:dim(SET_1_spectra_D2_SPLINE)[2]]
colnames(SET_2_spectra_D2_SPLINE) = row_names_all[1:dim(SET_2_spectra_D2_SPLINE)[2]]
rownames(rf$importance) = row_names_all[1:length(rf$importance)]
varImpPlot(rf, n.var = 15)

# CORRELATION BETWEEN WAVENUMBERS AND OUTPUT
require(corrplot)
fiveteen_variables=order(rf$importance, decreasing = TRUE)[1:10]
##New data for correlogramme
Imp_var_data = cbind(SET_1_spectra_D2_SPLINE[,fiveteen_variables],data.frame(INFLAM_SET1$NRF2[-na],INFLAM_SET1$A2M[-na],INFLAM_SET1$Crp[-na]))
colnames(Imp_var_data) = c(colnames(SET_1_spectra_D2_SPLINE[,fiveteen_variables]),'Nrf2','A2m', 'Crp')
Imp_var_data_val = cbind(SET_2_spectra_D2_SPLINE[,fiveteen_variables],INFLAM_SET2)
colnames(Imp_var_data_val) = c(colnames(SET_2_spectra_D2_SPLINE[,fiveteen_variables]),'Nrf2','A2m', 'Crp')

##CORRELOGRAMME
r=rbind(Imp_var_data, Imp_var_data_val)
C.Imp = data.frame(cor(rbind(Imp_var_data, Imp_var_data_val)))
new = as.matrix(C.Imp[c(1:10), c(11:13)])
rownames(new) = c(colnames(Imp_var_data)[1:10])
colnames(new) = c('Nrf2', "A2m", 'Crp')
##PLOT OF CORRELOGRAMME
library(plot.matrix)
library(RColorBrewer)
par(fig=c(0,0.6,0.02,0.98))
plot(new, main='', xlab='', ylab='', yaxt='n',xaxt = "n" ,cex.main=0.7,
     border=NA, text.cell=list(cex=0.7), digits=2, 
     col=brewer.pal(n = 9, name = 'BuGn'), key=NULL)
##PLOT OF R SQUARED OF EACH WAVENUMBERS AND A2M.
R_squared=R_Squared2=p_value=NULL
for(i in 1:10){
  R=summary(lm(r$A2m~r[,i]))
  R_squared=c(R_squared,R$r.squared)
  p_value=c(p_value,R$coefficients[2,4])
}
par(fig=c(0.40,1,0,1), new=TRUE)
barplot(rev(R_squared),col="grey", horiz = TRUE, main='', density=20)
