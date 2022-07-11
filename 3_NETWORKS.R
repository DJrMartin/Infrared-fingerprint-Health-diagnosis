installed = installed.packages()
if (!("igraph" %in% installed)) { install.packages("igraph") }

library('igraph')
#DATA IMPORT FROM RANDOM FOREST IMPORTANT WAVENUMBERS

# GRAPH FOR STEATOSIS/TRIGLY
graph='steatosis'
if (graph=='trigly'){
  load('~/Dropbox/DPSMET_2021/codes figures papier/networks_for_trigly.RData')
}else{
  load('~/Dropbox/DPSMET_2021/codes figures papier/networks_for_steatosis.RData')
}

net <- graph_from_data_frame(d=links, vertices=nodes, directed=T) 
net

colrs <- c("gold", "tomato", "lightgreen", "green3")
V(net)$color <- colrs[V(net)$type]
V(net)$size <- abs(log(V(net)$importance+1))*39
plot(net, edge.arrow.size=0, vertex.label.cex=0.8,vertex.label.color="black",margin=0,
     mark.col=c("#C5E5E7","#ECD89A"), mark.border=NA)
legend(x=-1.5, y=-1.1, c("C=O (Esters)", "C-O (Lipids)","C-O (Carbohydrate)","C-C/C-O (Carbohydrate)"), pch=21,
       col=colrs, pt.bg=colrs, pt.cex=2, cex=.8, bty="n", ncol=1)

# GRAPH FOR INFLAMMATION
load('~/Dropbox/DPSMET_2021/codes figures papier/network_for_inflammation.RData')

net <- graph_from_data_frame(d=links, vertices=nodes, directed=T) 
net
colrs <- c("gold","lightgreen","yellowgreen", "green3")

V(net)$color <- colrs[V(net)$type]
V(net)$size <- abs(log(V(net)$importance+1))*39

plot(net, edge.arrow.size=0, vertex.label.cex=0.8,vertex.label.color="black",margin=0,
     mark.col=c("#C5E5E7","#ECD89A"), mark.border=NA)

legend(x=-1.5, y=-1.1, c("C=O (Esters)","Amide I","Amide II", "Amide III"), pch=21,
       col=colrs, pt.bg=colrs, pt.cex=2, cex=.8, bty="n", ncol=1)
