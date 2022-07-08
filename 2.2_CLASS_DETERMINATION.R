setwd(dir="~/Desktop/DATA_HEALTH_DIAGNOSIS")
dir()
MIR.1 = read.csv("Blood_smear_bruts_SET1.csv", sep=';' )
MIR.2 = read.csv("Blood_smear_bruts_SET2.csv", sep=';')

#DATA OF INFLAMMATION MARKERS IN LIVER
INFLAM_SET1 = data.frame(scale(MIR.1[, c(19,20,22)]))
INFLAM_SET2 = data.frame(scale(MIR.2[, c(14:16)]))

na=which(is.na(INFLAM_SET1$NRF2))

#CLASSIFICATION OF LIVER INFLAMMATION FROM NRF2, CRP and A2M
md=dist(INFLAM_SET1[-na,])
arbre <- hclust(md, method = "ward.D2")
plot(arbre)

##DETERMINATION OF INFLAMMATION STATE IN THE SET 2
group_set2= (INFLAM_SET2$A2M.1>-0.18)&(INFLAM_SET2$CRP.1>-0.3)
INFLAM_SET2$Inflammation=group_set2
INFLAM_SET2$Inflammation[group_set2==TRUE]="Moderate"
INFLAM_SET2$Inflammation[group_set2==FALSE]="Low"
INFLAM_SET2$Inflammation=as.factor(INFLAM_SET2$Inflammation)

#VISUALISATION OF POSITION OF EACH MICE IN FUNCTION OF CRP AND A2M.
group = cutree(arbre,2)
w=which(group_set2==TRUE)
plot(INFLAM_SET1[-na,2:3],col=group,pch=20,
     xlim = range(INFLAM_SET1[-na,2],INFLAM_SET2[,2]),
     ylim = range(INFLAM_SET1[-na,3],INFLAM_SET2[,3]))
points(INFLAM_SET2[w,2:3],pch=2,col="blue")
points(INFLAM_SET2[-w,2:3],pch=2,col="Green3")
legend('topleft', legend=c('Low', "Moderate"), fill=c('red', 'black'), bty='n', title = 'SET 1', cex=0.7)
legend('bottomright', legend=c('Low', "Moderate"), fill=c('Green3', 'blue'), bty='n', title = 'SET 2', cex=0.7)

###Figure S5.1
palette <- hcl.colors(50, "Spectral", rev = TRUE)
heatmap(as.matrix(INFLAM_SET1[-na,]), Colv=T,hclustfun=function(x) hclust(x, method="ward.D2"), scale='none', col=palette,cexCol=0.8)
legend(x="topleft", legend=c("min","average","max"), 
       fill=hcl.colors(50, "Spectral", rev = TRUE)[c(1,25,50)])

##Figure S5.2
layout(matrix(c(1,1,
                2,2), nrow=2, byrow=TRUE))
pie(table(MIR.1$Group[which(group==2)]), col=c(hcl.colors(4, "Spectral", rev = TRUE)), main='LOW \n INFLAMMATION')
pie(table(MIR.1$Group[which(group==1)]), col=c(hcl.colors(4, "Spectral", rev = TRUE)), main="MODERATE \n INFLAMMATION")
