library(signal)
library(pls)
library(fda)
library(FactoMineR)

#DATA IMPORTATION
setwd(dir="~/Desktop/DATA_HEALTH_DIAGNOSIS")
dir()
MIR.1 = read.csv("Serum_bruts_SET1.csv", sep=';')
MIR.2 = read.csv("Serum_bruts_SET2.csv", sep=';')

#WAVENUMBERS SELECTION
freq_set1=as.numeric(substr(colnames(MIR.1),2,9))
freq_set2=as.numeric(substr(colnames(MIR.2),2,9))

fingerprint_set1=c(which(freq_set1<1800&freq_set1>800))
fingerprint_set2=c(which(freq_set2<1800&freq_set2>800))

#DATA VISUALISATION
matplot(freq_set1[fingerprint_set1], t(MIR.1[,fingerprint_set1]), typ="l",col=1,lty=1,ylab="Absorbances",xlab="freq")
matlines(freq_set2[fingerprint_set2], t(MIR.2[,fingerprint_set2]), typ="l",col=2,lty=1)

SET_1_spectra=MIR.1[,fingerprint_set1]
SET_2_spectra=MIR.2[,fingerprint_set2]

dev.off()

#DATA TRANSFORMATIONS, SECOND DERIVES
##SET 1
Sd2raman = apply(SET_1_spectra,1,sgolayfilt,p=3,n=9,m=2)
str(Sd2raman)
Sd2raman = t(Sd2raman)
colnames(Sd2raman)=freq_set1[fingerprint_set1]
SET_1_spectra_D2=as.data.frame(Sd2raman)
##SET 2
Sd2raman = apply(SET_2_spectra,1,sgolayfilt,p=3,n=9,m=2)
str(Sd2raman)
Sd2raman = t(Sd2raman)
colnames(Sd2raman)=freq_set1[fingerprint_set1]
SET_2_spectra_D2=as.data.frame(Sd2raman)

matplot(freq_set1[fingerprint_set1], t(SET_1_spectra_D2), typ="l",col=1,lty=1,ylab="Absorbances",xlab="freq")
matlines(freq_set1[fingerprint_set1], t(SET_2_spectra_D2), typ="l",col=2,lty=1)

dev.off()

#DATA NORMALISATION THROUGH VECTOR NORMALISATION
##Normalise vector
normalize.vector <-function(v){v/sqrt(sum(v^2))}
MAT.D2=rbind(SET_1_spectra_D2,SET_2_spectra_D2)
all_spec <- normalize.vector(MAT.D2)

SET_1_spectra_D2=as.data.frame(all_spec[c(1:93),])
SET_2_spectra_D2=as.data.frame(all_spec[c(94:119),])

#OUTLIERS IDENTIFICATION
SET_1_2_SPECTRA=data.frame(rbind(SET_1_spectra_D2,SET_2_spectra_D2))
res.PCA=PCA(SET_1_2_SPECTRA)
ylim=range(res.PCA$ind$coord[,2])
xlim=range(res.PCA$ind$coord[,1])
plot(data.frame(res.PCA$ind$coord[c(1:93),c(1,2)]), col=1, ylim=ylim, xlim=xlim,
     xlab='PC1', ylab='PC2')
points(res.PCA$ind$coord[c(94:119),c(1,2)], col=2)

paste('individu',c(which(as.numeric(res.PCA$ind$coord[,1])>40),
                   which(as.numeric(res.PCA$ind$coord[,1])<(-40))),
      'from set 1')
#33, 35,49, 88
outliers=c(33,35,49,88)

#DATA COMPRESSION THROUGH A B SPLINE FUNCTION
projRecomp = function(xdata, nBasis, t = 1:dim(xdata)[2], basis = "Splines"){
  t = sort((t-min(t))/(max(t)-min(t)))
  if (basis == "Fourier") {basisobj = create.fourier.basis(nbasis = nBasis)}
  if (basis == "Splines") {basisobj = create.bspline.basis(norder = 3, breaks = seq(head(t,1), tail(t,1), length = nBasis-1))}
  BFunction = getbasismatrix(t, basisobj) 
  Fdata = t(sapply(1:nrow(xdata), 
                   FUN = function(i) t(solve(t(BFunction)
                                             %*% BFunction)%*% t(BFunction) %*% xdata[i,])))
  FdataRec = t(sapply(1:nrow(xdata), 
                      FUN = function(i) t(BFunction %*% solve(t(BFunction)%*% BFunction)%*% 
                                            t(BFunction) %*% xdata[i,])))
  return(list(coeffProj = Fdata, foncRecon = FdataRec, BFunction = BFunction, basisobj = basisobj))
}

b=c(7,5,2)
p = 475
VI = matrix(0,length(b),p)
x_all = NULL
cnt = 1
L = 0
x_recom=NULL
for (i in b){
  nBasis = round(p/i)
  L = L + nBasis
  x = projRecomp(as.matrix(SET_1_2_SPECTRA), nBasis = nBasis)
  print(nBasis)
  print(dim(x$coeffProj))
  if (cnt==1){
    x_all = x$coeffProj
    x_recom= x$foncRecon
    print(dim(x_recom))
  } else {
    x_all = cbind(x_all,x$coeffProj)
    x_recom= cbind(x_recom, x$foncRecon)
    print(dim(x_all))
    print(dim(x_recom))
  }
  cnt = cnt+1
}
plot(freq_set1[fingerprint_set1][1:485],SET_1_2_SPECTRA[1,],typ="l",
     xlab="frequence (in cm-1)",ylab="Absorbance", ylim=range(x_all[1,]))
  lines(freq_set1[fingerprint_set1][1:485],x_recom[1,1:485],col="red")
  lines(freq_set1[fingerprint_set1][1:485],x_recom[1,486:970],col="green")
  lines(freq_set1[fingerprint_set1][1:485],x_recom[1,971:1455],col="blue")

#DATA SAVING FOR NEXT EXPLORATION
SET_1_spectra_D2_SPLINE=as.data.frame(x_all[c(1:93),])
SET_1_spectra_D2_SPLINE=as.data.frame(SET_1_spectra_D2_SPLINE[-outliers,])
SET_2_spectra_D2_SPLINE=as.data.frame(x_all[c(94:119),])

save(SET_1_spectra_D2_SPLINE,
     SET_2_spectra_D2_SPLINE,
     file='SET_1_2.RData')
