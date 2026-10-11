#!/usr/bin/env Rscript
# Cluster UMAPs and their marker rows share one vertical block layout.
args0<-commandArgs(FALSE);self<-sub('^--file=','',args0[grepl('^--file=',args0)])
self<-gsub('~+~',' ',self,fixed=TRUE);ROOT<-normalizePath(file.path(dirname(self),'..'))
args<-commandArgs(TRUE);if(length(args)<2)stop('Usage: render_marker_atlas.R marker_atlas.rds output_dir [prefix]')
input<-normalizePath(args[1]);OUT<-args[2];prefix<-if(length(args)>2)args[3] else 'cluster_markers'
for(p in c('grid','ragg','svglite','systemfonts'))if(!requireNamespace(p,quietly=TRUE))stop('Missing package: ',p)
source(file.path(ROOT,'assets/style/approved_style.R'));library(grid)
b<-readRDS(input);cells<-b$cells;cl<-b$clusters;mk<-b$markers
if(!'patient_id'%in%names(cells))cells$patient_id<-''
title<-if(is.null(b$title))'Cluster markers' else b$title
stopifnot(!anyDuplicated(cells$cell_id),!anyDuplicated(mk$row_id),
          identical(as.character(mk$feature),rownames(b$dot_mean_z)),
          identical(rownames(b$cell_z),as.character(mk$feature)),ncol(b$cell_z)==nrow(cells),
          identical(colnames(b$cell_z),cells$cell_id),setequal(cl,cells$cluster_id),
          all(is.finite(b$cell_z)),all(is.finite(b$dot_mean_z)),all(b$dot_pct>=0&b$dot_pct<=100))
if(!STYLE$font%in%systemfonts::system_fonts()$family)stop('Arial unavailable')
dir.create(OUT,recursive=TRUE,showWarnings=FALSE)
mm<-function(x)unit(x,'mm')
txt<-function(s,x,y,pt=8,just='centre',rot=0,col='#202020')grid.text(s,mm(x),mm(y),just=just,rot=rot,gp=gpar(fontfamily=STYLE$font,fontsize=pt,col=col))
rect<-function(x,y,w,h,col)grid.rect(mm(x),mm(y),mm(w),mm(h),gp=gpar(fill=col,col=NA))
circle<-function(x,y,d,col)grid.circle(mm(x),mm(y),r=mm(d/2),gp=gpar(fill=col,col=NA))
vp<-function(x,y,w,h)pushViewport(viewport(x=mm(x),y=mm(y),width=mm(w),height=mm(h),just=c('left','bottom'),clip='off'))
rh<-5.7;gap<-7;gh<-pmax(34,as.integer(table(factor(mk$group,levels=cl)))*rh)
bottom<-29;top<-bottom+sum(gh)+gap*(length(cl)-1);centres<-ys<-numeric()
cursor<-top;row_y<-numeric(nrow(mk))
for(k in cl){ix<-which(mk$group==k);h<-gh[match(k,cl)];centres[k]<-cursor-h/2;row_y[ix]<-centres[k]+((length(ix)+1)/2-seq_along(ix))*rh;cursor<-cursor-h-gap}
mini_x<-24;gene_x<-66;dot_x<-seq(77,105,length.out=length(cl));heat_start<-119;heat_w<-23;heat_gap<-1
W<-heat_start+length(cl)*(heat_w+heat_gap)+8;hero_size<-135;hero_y<-top+34;H<-hero_y+hero_size+17
span<-max(diff(range(cells$umap_1)),diff(range(cells$umap_2)))*1.10
xmin<-mean(range(cells$umap_1))-span/2;ymin<-mean(range(cells$umap_2))-span/2
ux<-function(x)(x-xmin)/span;uy<-function(y)(y-ymin)/span
set.seed(194247);draw_order<-sample(seq_len(nrow(cells)));factor<-max(.35,min(1.5,sqrt(11941/nrow(cells))))
cloud<-function(k=NULL){
 cap<-ragg::agg_capture(width=if(is.null(k))110 else 30,height=if(is.null(k))110 else 30,units='mm',res=600,background='white');grid.newpage()
 if(!is.null(k))grid.points(ux(cells$umap_1),uy(cells$umap_2),default.units='npc',pch=16,size=mm(.20*factor),gp=gpar(col=adjustcolor('#C8C8C8',alpha.f=.72)))
 ix<-if(is.null(k))draw_order else draw_order[cells$cluster_id[draw_order]==k]
 grid.points(ux(cells$umap_1[ix]),uy(cells$umap_2[ix]),default.units='npc',pch=16,size=mm(if(is.null(k)).43*factor else .20*factor),gp=gpar(col=adjustcolor(b$colors[cells$cluster_id[ix]],alpha.f=.85)))
 im<-cap();dev.off();im
}
hero<-cloud();minis<-setNames(lapply(cl,cloud),cl)
labs<-do.call(rbind,lapply(cl,function(k){d<-cells[cells$cluster_id==k,];i<-which.min((d$umap_1-median(d$umap_1))^2+(d$umap_2-median(d$umap_2))^2);data.frame(cluster_id=k,x=d$umap_1[i],y=d$umap_2[i])}))
ramps<-read.csv(file.path(ROOT,'assets/style/continuous_scales.csv'))
dot_ramp<-ramps$color[ramps$scale=='or'];heat_ramp<-ramps$color[ramps$scale=='pathway']
map<-function(x,ramp,limit){t<-pmax(0,pmin(1,(x+limit)/(2*limit)));ramp[round(t*(length(ramp)-1))+1]}
# Every cluster's original cells are included, sorted by donor and barcode.
cell_orders<-setNames(lapply(cl,function(k){ix<-which(cells$cluster_id==k);ix[order(cells$patient_id[ix],cells$cell_id[ix])]}),cl)
heat_images<-list()
for(i in seq_along(cl))for(j in seq_along(cl)){
 ix<-which(mk$group==cl[i]);cc<-cell_orders[[cl[j]]]
 mat<-matrix(map(b$cell_z[ix,cc,drop=FALSE],heat_ramp,2),nrow=length(ix),ncol=length(cc))
 heat_images[[paste(i,j)]]<-as.raster(mat)
}
bar<-function(ramp,x,y,title,limit){
 w<-22;for(i in seq_along(ramp))rect(x-w/2+(i-.5)*w/length(ramp),y,w/length(ramp)+.01,2,ramp[i])
 txt(title,x,y+4,7);txt(paste0('<=-',limit),x-11,y-4,6.5);txt('0',x,y-4,6.5);txt(paste0('>=',limit),x+11,y-4,6.5)
}
draw<-function(){
 grid.newpage();txt(title,W/2,H-6,10)
 vp((W-hero_size)/2,hero_y,hero_size,hero_size);grid.raster(hero,width=unit(1,'npc'),height=unit(1,'npc'),interpolate=FALSE)
 for(i in seq_len(nrow(labs))){xx<-ux(labs$x[i])*hero_size;yy<-uy(labs$y[i])*hero_size
  for(a in seq(0,2*pi,length.out=13)[-13])txt(labs$cluster_id[i],xx+.19*cos(a),yy+.19*sin(a),8.6,col='white')
  txt(labs$cluster_id[i],xx,yy,8.6)
 }
 grid.lines(mm(c(3,3,18)),mm(c(18,4,4)),gp=gpar(lwd=.55));txt('UMAP-1',11,1,7);txt('UMAP-2',.5,11,7,rot=90)
 txt(paste0('n = ',format(nrow(cells),big.mark=',')),23,6,7,just='left');popViewport()
 txt('Cluster',mini_x,top+17,8);txt('Markers',gene_x-6,top+17,8)
 txt('Expression by cluster',mean(dot_x),top+17,8)
 txt('Single-cell expression',heat_start+(length(cl)*(heat_w+heat_gap)-heat_gap)/2,top+17,8)
 for(j in seq_along(cl)){
  txt(cl[j],dot_x[j],top+5,7.5);circle(dot_x[j],top+9,1.5,b$colors[cl[j]])
  x<-heat_start+(j-1)*(heat_w+heat_gap)+heat_w/2
  txt(cl[j],x,top+7,7.5);txt(format(length(cell_orders[[j]]),big.mark=','),x,top+3,6,col='#666666');rect(x,top+.7,heat_w,1,b$colors[cl[j]])
 }
 for(i in seq_along(cl)){
  k<-cl[i];ix<-which(mk$group==k);cy<-centres[k]
  grid.raster(minis[[k]],mm(mini_x),mm(cy),width=mm(29),height=mm(29),interpolate=FALSE)
  txt(k,mini_x,cy+16.8,8.5)
  for(r in ix){
   txt(mk$feature[r],gene_x,row_y[r],8,just='right')
   for(j in seq_along(cl))circle(dot_x[j],row_y[r],3.5*sqrt(b$dot_pct[r,j]/100),map(b$dot_mean_z[r,j],dot_ramp,1.5))
  }
  for(j in seq_along(cl)){
   x<-heat_start+(j-1)*(heat_w+heat_gap)+heat_w/2
   grid.raster(heat_images[[paste(i,j)]],mm(x),mm(cy),width=mm(heat_w),height=mm(length(ix)*rh),interpolate=FALSE)
  }
 }
 bar(dot_ramp,mean(dot_x),17,'Mean expression (z)',1.5)
 txt('Detected cells (%)',52,21,7)
 for(j in 1:3){v<-c(25,50,100)[j];xx<-43+(j-1)*9;circle(xx,16,3.5*sqrt(v/100),'#BDBDBD');txt(as.character(v),xx,11,6.5)}
 bar(heat_ramp,heat_start+(length(cl)*(heat_w+heat_gap)-heat_gap)/2,17,'Cell expression (z)',2)
 txt('Heatmap: all cells; equal-width cluster blocks',heat_start+(length(cl)*(heat_w+heat_gap)-heat_gap)/2,4,6,col='#666666')
}
stem<-file.path(OUT,prefix)
ragg::agg_png(paste0(stem,'.png'),width=W,height=H,units='mm',res=450);draw();dev.off()
if(Sys.info()[['sysname']]=='Darwin')quartz(type='pdf',file=paste0(stem,'.pdf'),width=W/25.4,height=H/25.4,family=STYLE$font) else cairo_pdf(paste0(stem,'.pdf'),width=W/25.4,height=H/25.4,family=STYLE$font)
draw();dev.off();svglite::svglite(paste0(stem,'.svg'),width=W/25.4,height=H/25.4);draw();dev.off()
write.table(data.frame(row_id=mk$row_id,cluster=mk$group,gene=mk$feature,y_mm=row_y),file.path(OUT,'marker_row_alignment.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
write.table(data.frame(cluster=cl,n_cells=lengths(cell_orders),mini_y_mm=centres[cl]),file.path(OUT,'highlight_alignment.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
capture.output(dput(list(input_md5=tools::md5sum(input),metadata=b$metadata,width_mm=W,height_mm=H)),file=file.path(OUT,'render_record.R'))
cat('Rendered marker-first atlas:',nrow(mk),'marker rows,',nrow(cells),'cells.\n')
