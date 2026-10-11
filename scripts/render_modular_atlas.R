#!/usr/bin/env Rscript
# Accepted vertical scaffold. Biology and optional modules come only from the job RDS.
suppressPackageStartupMessages({library(grid);library(ragg);library(svglite)})
a<-commandArgs(TRUE);if(length(a)<3)stop('Usage: render_modular_atlas.R figure_job.rds output_dir stem')
J<-readRDS(a[1]);O<-a[2];dir.create(O,recursive=TRUE,showWarnings=FALSE);stem<-a[3]
self<-sub('^--file=','',grep('^--file=',commandArgs(FALSE),value=TRUE)[1]);root<-dirname(dirname(normalizePath(self)))
m<-J$marker;cl<-as.character(m$clusters);K<-length(cl);ps<-J$panels;if(is.null(ps))ps<-list();S<-J$style
setting<-function(n,d)if(is.null(S[[n]]))d else S[[n]]
W<-setting('width_mm',if(length(ps))183 else 100);hs<-setting('hero_mm',min(98,W-10));band<-setting('band_mm',17.5);zlim<-setting('zlim',2);pd<-setting('marker_point_mm',1.8)
stopifnot(K>=2,K<=9,nrow(m$markers)==nrow(m$dot_mean_z),nrow(m$markers)==nrow(m$dot_pct),identical(colnames(m$dot_mean_z),cl),identical(colnames(m$dot_pct),cl),all(is.finite(m$dot_mean_z)),all(is.finite(m$dot_pct)),all(m$dot_pct>=0&m$dot_pct<=100),all(m$cells$cluster_id%in%cl),all(table(m$markers$group)<=6))
cy<-25.5+(rev(seq_len(K))-1)*band;title_y<-cy[1]+28;hero_bottom<-title_y+4;hy<-hero_bottom+hs/2;H<-hero_bottom+hs+4;hx<-W/2+2.5
if(H>247)message('Figure height is ',H,' mm: review target journal dimensions before submission.')
mm<-function(x)unit(x,'mm')
txt<-function(s,x,y,pt=6.5,just='centre',rot=0,col='#222222',face='plain')grid.text(s,mm(x),mm(y),just=just,rot=rot,gp=gpar(fontfamily='Arial',fontsize=pt,col=col,fontface=face))
rect<-function(x,y,w,h,c)grid.rect(mm(x),mm(y),mm(w),mm(h),gp=gpar(fill=c,col=NA))
circ<-function(x,y,d,c)grid.circle(mm(x),mm(y),r=mm(d/2),gp=gpar(fill=c,col=NA))
ln<-function(x,y,col='#555555',lwd=.5)grid.lines(mm(x),mm(y),gp=gpar(col=col,lwd=lwd))
tag<-function(s,x,y)txt(s,x,y,8,face='bold',just='left')
r<-read.csv(file.path(root,'assets/style/measured_Fig4_5_colorbars.csv'));ramp<-r$hex[r$figure==4&r$scale=='pathway']
colz<-function(z)ramp[1+round(pmax(0,pmin(1,(z+zlim)/(2*zlim)))*(length(ramp)-1))]
span<-max(diff(range(m$cells$umap_1)),diff(range(m$cells$umap_2)))*1.1
ux<-function(x)(x-mean(range(m$cells$umap_1))+span/2)/span;uy<-function(y)(y-mean(range(m$cells$umap_2))+span/2)/span
set.seed(setting('seed',1));ord<-sample(seq_len(nrow(m$cells)))
cloud<-function(k=NULL,sz=hs){cap<-ragg::agg_capture(width=sz,height=sz,units='mm',res=900,background='white');grid.newpage()
 if(!is.null(k))grid.points(ux(m$cells$umap_1),uy(m$cells$umap_2),default.units='npc',pch=16,size=mm(.105),gp=gpar(col=adjustcolor('#C9C9C9',alpha.f=.7)))
 ii<-if(is.null(k))ord else ord[m$cells$cluster_id[ord]==k]
 grid.points(ux(m$cells$umap_1[ii]),uy(m$cells$umap_2[ii]),default.units='npc',pch=16,size=mm(if(is.null(k)).36 else .12),gp=gpar(col=adjustcolor(m$colors[m$cells$cluster_id[ii]],alpha.f=.85)))
 z<-cap();dev.off();z}
hero<-cloud();minis<-setNames(lapply(cl,function(k)cloud(k,15)),cl)
# Prepare geometry and validate module data before drawing.
left<-80
for(q in seq_along(ps)){p<-ps[[q]];stopifnot(p$type%in%c('heatmap','violin','replicate_points'),length(p$features)>0,!anyDuplicated(p$features),p$width_mm>0)
 p$left<-left;p$cell_width<-p$width_mm/length(p$features);p$x<-left+(seq_along(p$features)-.5)*p$cell_width
 if(!is.null(p$groups)){
  stopifnot(identical(as.integer(unlist(p$groups,use.names=FALSE)),seq_along(p$features)))
  gap<-if(is.null(p$group_gap_mm))2.5 else p$group_gap_mm
  p$cell_width<-(p$width_mm-(length(p$groups)-1)*gap)/length(p$features);stopifnot(p$cell_width>0)
  gi<-rep(seq_along(p$groups),lengths(p$groups));p$x<-left+(seq_along(p$features)-.5)*p$cell_width+(gi-1)*gap
 }
 if(p$type%in%c('heatmap','violin'))stopifnot(identical(rownames(p$values),p$features),identical(colnames(p$values),cl),all(is.finite(p$values)))
 if(p$type=='violin'){
  stopifnot(identical(rownames(p$expression),p$features),identical(colnames(p$expression),m$cells$cell_id),all(is.finite(p$expression)),all(p$expression>=0))
  if(is.null(p$ylim))p$ylim<-c(0,ceiling(max(p$expression)*2)/2)
  stopifnot(p$ylim[1]==0,p$ylim[2]>=max(p$expression));p$density<-list()
  for(g in p$features){bw<-max(.12,bw.nrd0(p$expression[g,]));for(k in cl){v<-p$expression[g,m$cells$cluster_id==k];d<-density(c(v,-v),bw=bw,from=0,to=p$ylim[2],n=256,cut=0);p$density[[paste(g,k)]]<-list(y=d$x,w=d$y/max(d$y))}}
 }
 if(p$type=='replicate_points'){
  stopifnot(all(c('cluster','replicate_id','gene','mean')%in%names(p$data)),all(p$data$cluster%in%cl),all(p$data$gene%in%p$features),all(is.finite(p$data$mean)),!anyDuplicated(paste(p$data$cluster,p$data$replicate_id,p$data$gene)))
  stopifnot(length(p$ylim)==2,p$ylim[1]<p$ylim[2],all(p$data$mean>=p$ylim[1]&p$data$mean<=p$ylim[2]))
  p$data<-p$data[order(p$data$cluster,p$data$gene,p$data$replicate_id),];set.seed(setting('seed',1)+q);p$data$jitter_mm<-runif(nrow(p$data),-2.7,2.7)
 }
 ps[[q]]<-p;left<-left+p$width_mm+7
}
if(length(ps)&&left-7>W-5)stop('Configured modules exceed canvas width; revise panel widths or canvas.')
marker_x<-seq(51,71.5,length.out=K)
draw<-function(){grid.newpage();tag('a',3,H-5)
 grid.raster(hero,mm(hx),mm(hy),width=mm(hs),height=mm(hs),interpolate=FALSE)
 for(k in cl){d<-m$cells[m$cells$cluster_id==k,];j<-which.min((d$umap_1-median(d$umap_1))^2+(d$umap_2-median(d$umap_2))^2);x<-hx-hs/2+ux(d$umap_1[j])*hs;y<-hy-hs/2+uy(d$umap_2[j])*hs
  for(th in seq(0,2*pi,length.out=13)[-13])txt(k,x+.13*cos(th),y+.13*sin(th),7,col='white');txt(k,x,y,7)}
 ax<-hx-hs/2+2;ay<-hero_bottom+4
 grid.lines(mm(c(ax,ax)),mm(c(ay,ay+12)),arrow=arrow(length=mm(1.1),type='open'),gp=gpar(col='#555555',lwd=.5))
 grid.lines(mm(c(ax,ax+12)),mm(c(ay,ay)),arrow=arrow(length=mm(1.1),type='open'),gp=gpar(col='#555555',lwd=.5))
 txt('UMAP 1',ax+6,ay-3,6);txt('UMAP 2',ax-3.5,ay+6,6,rot=90)
 tag('b',3,title_y);txt('Markers',mean(marker_x),title_y,7)
 for(j in seq_along(cl)){txt(cl[j],marker_x[j],cy[1]+16,6);circ(marker_x[j],cy[1]+19,1.2,m$colors[cl[j]])}
 for(p in ps){
  if(is.null(p$groups))txt(p$title,mean(p$x),title_y,6.5) else for(nm in names(p$groups)){ii<-p$groups[[nm]];txt(nm,mean(p$x[ii]),title_y,6.5);ln(range(p$x[ii])+c(-p$cell_width/2,p$cell_width/2),rep(title_y-3,2))}
  for(j in seq_along(p$features))txt(p$features[j],p$x[j],cy[1]+11,6.5,just='left',rot=55)
 }
 for(i in seq_along(cl)){k<-cl[i];y<-cy[i];grid.raster(minis[[k]],mm(14),mm(y),width=mm(15),height=mm(15),interpolate=FALSE);txt(k,14,y+8.6,6.5)
  ii<-which(m$markers$group==k)
  for(j in seq_along(ii)){yy<-y+((length(ii)+1)/2-j)*2.6;txt(m$markers$feature[ii[j]],46,yy,6,just='right');for(z in seq_along(cl))circ(marker_x[z],yy,pd*sqrt(m$dot_pct[ii[j],z]/100),colz(m$dot_mean_z[ii[j],z]))}
  for(p in ps)for(j in seq_along(p$features)){g<-p$features[j];x<-p$x[j]
   if(p$type=='heatmap')rect(x,y,p$cell_width+.01,7.4,colz(p$values[g,k]))
   if(p$type=='violin'){d<-p$density[[paste(g,k)]];yy<-y-4.9+d$y/p$ylim[2]*9.8;ww<-d$w*min(3,p$width_mm/length(p$features)*.36);grid.polygon(mm(c(x-ww,rev(x+ww))),mm(c(yy,rev(yy))),gp=gpar(fill=colz(p$values[g,k]),col=NA))}
   if(p$type=='replicate_points'){
    d<-p$data[p$data$cluster==k&p$data$gene==g,];ln(c(x-3.5,x+3.5),c(y-4.9,y-4.9),col='#DDDDDD',lwd=.35)
    if(nrow(d)){yy<-y-4.9+(d$mean-p$ylim[1])/diff(p$ylim)*9.8;for(h in seq_len(nrow(d)))circ(x+d$jitter_mm[h],yy[h],.65,adjustcolor(m$colors[k],alpha.f=.72));ym<-y-4.9+(median(d$mean)-p$ylim[1])/diff(p$ylim)*9.8;ln(c(x-2.8,x+2.8),rep(ym,2),col='#333333',lwd=.65)}
   }
  }
 }
 txt('Detected (%)',62,12,6);for(j in 1:3){v<-c(25,50,100)[j];x<-55+(j-1)*7;circ(x,7.8,pd*sqrt(v/100),'#B7B7B7');txt(v,x,4,6)}
 lx<-if(length(ps))102 else 28;txt('Mean expression (z)',lx,12,6)
 for(j in seq_along(ramp))rect(lx-11+(j-.5)*22/length(ramp),7.8,22/length(ramp)+.035,1.7,ramp[j])
 txt(paste0('≤−',zlim),lx-11,4,6);txt('0',lx,4,6);txt(paste0('≥',zlim),lx+11,4,6)
 for(p in ps)if(p$type=='replicate_points'){
  x<-p$left+p$width_mm+1;y<-cy[1];ln(c(x,x),c(y-4.9,y+4.9),lwd=.4)
  for(v in p$ylim){yy<-y-4.9+(v-p$ylim[1])/diff(p$ylim)*9.8;ln(c(x,x+1),c(yy,yy),lwd=.4);txt(format(v),x+1.8,yy,6,just='left')}
  txt(p$unit,mean(p$x),12,6);circ(mean(p$x)-7,6, .65,'#666666');txt(if(is.null(p$replicate_label))'Sample' else p$replicate_label,mean(p$x)-5.5,6,6,just='left')
 }
}
f<-file.path(O,stem)
ragg::agg_png(paste0(f,'.png'),width=W,height=H,units='mm',res=600,background='white');draw();dev.off()
if(capabilities('aqua'))grDevices::quartz(type='pdf',file=paste0(f,'.pdf'),width=W/25.4,height=H/25.4,family='Arial') else grDevices::cairo_pdf(paste0(f,'.pdf'),width=W/25.4,height=H/25.4,family='Arial')
draw();dev.off();svglite::svglite(paste0(f,'.svg'),width=W/25.4,height=H/25.4);draw();dev.off()
ragg::agg_png(paste0(f,'_preview.png'),width=W,height=H,units='mm',res=220,background='white');draw();dev.off()
write.table(data.frame(cluster=cl,y_mm=cy,marker_n=as.integer(table(factor(m$markers$group,levels=cl)))),file.path(O,paste0(stem,'_alignment.tsv')),sep='\t',quote=FALSE,row.names=FALSE)
saveRDS(list(width_mm=W,height_mm=H,clusters=cl,panels=lapply(ps,function(p)p[c('type','features','ylim','x','title')]),zlim=zlim),file.path(O,paste0(stem,'_render_metadata.rds')))
cat('Rendered',K,'clusters;',nrow(m$markers),'marker rows;',length(ps),'optional modules;',W,'x',H,'mm.\n')
