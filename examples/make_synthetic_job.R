#!/usr/bin/env Rscript
# Synthetic coordinates and artificial features for rendering only.
set.seed(20261011)
args<-commandArgs(FALSE);self<-sub('^--file=','',grep('^--file=',args,value=TRUE)[1]);dest<-file.path(dirname(normalizePath(self)),'generated');dir.create(dest,showWarnings=FALSE)
cl<-paste0('C',0:5);n<-rep(1500,6);center<-rbind(c(0,-1.9),c(-.4,0),c(-1.7,1.1),c(1.5,2),c(.7,3),c(-1,3.4))
cells<-do.call(rbind,lapply(seq_along(cl),function(i){z<-matrix(rnorm(n[i]*2),ncol=2);theta<-i*.55;rot<-matrix(c(cos(theta),-sin(theta),sin(theta),cos(theta)),2);z<-sweep(z,2,c(.7,.32),'*')%*%rot;z<-sweep(z,2,center[i,],'+');data.frame(cell_id=paste0('synthetic_',i,'_',seq_len(n[i])),cluster_id=cl[i],umap_1=z[,1],umap_2=z[,2])}))
markers<-data.frame(group=rep(cl,each=4),feature=paste0('Feature_',seq_len(24)))
mean<-matrix(runif(24*6,.05,.5),24,6);pct<-matrix(runif(24*6,2,28),24,6)
for(i in seq_along(cl)){j<-which(markers$group==cl[i]);mean[j,i]<-runif(4,1.2,2);pct[j,i]<-runif(4,65,95)}
z<-t(scale(t(mean)));colnames(z)<-colnames(pct)<-cl
colors<-setNames(c('#EF7A84','#BE3A9D','#D0BE50','#80654B','#64BCC0','#7671CE'),cl)
heat<-matrix(c(1.2,.6,-.5,1.5,-.3,-1,.8,1.1,-.9,.1,.7,-1.5,-1,.2,1.4,.6,-.4,.9),3,6,dimnames=list(c('Axis_A','Axis_B','Axis_C'),cl))
expr<-matrix(rexp(nrow(cells)*2,rate=3),2,nrow(cells),dimnames=list(c('Signal_A','Signal_B'),cells$cell_id));expr[1,cells$cluster_id=='C3']<-expr[1,cells$cluster_id=='C3']+1;expr[2,cells$cluster_id=='C4']<-expr[2,cells$cluster_id=='C4']+1.4
v<-sapply(cl,function(k)rowMeans(expr[,cells$cluster_id==k,drop=FALSE]));v<-t(scale(t(v)));colnames(v)<-cl
job<-list(marker=list(cells=cells,clusters=cl,colors=colors,markers=markers,dot_mean_z=z,dot_pct=pct),panels=list(list(type='heatmap',title='Synthetic axes',features=rownames(heat),width_mm=39,values=heat),list(type='violin',title='Synthetic signals',features=rownames(expr),width_mm=43,values=v,expression=expr,ylim=c(0,ceiling(max(expr))))),style=list(width_mm=183))
saveRDS(job,file.path(dest,'figure_job.rds'))
writeLines('SYNTHETIC DATA: coordinates and all values are generated for a renderer demonstration; no biological inference.',file.path(dest,'DATA_NOTICE.txt'))
cat('Wrote synthetic job to',dest,'\n')
