#!/usr/bin/env Rscript
# Flexible rendering only. Input coordinates, assignments and values are unchanged.
args0<-commandArgs(FALSE);self<-sub("^--file=","",args0[grepl("^--file=",args0)])
self<-gsub("~+~"," ",self,fixed=TRUE)
ROOT<-normalizePath(file.path(dirname(self),".."),mustWork=TRUE)
args<-commandArgs(TRUE);if(!length(args))stop("Usage: Rscript render_atlas.R config.R [--validate-only]")
cfg<-normalizePath(args[1],mustWork=TRUE);e<-new.env(parent=baseenv());sys.source(cfg,e);A<-e$ATLAS
if(is.null(A))stop("Config must define ATLAS")
resolve<-function(p)if(grepl("^(/|[A-Za-z]:[/\\\\])",p))p else file.path(dirname(cfg),p)
DATA<-normalizePath(resolve(A$input_dir),mustWork=TRUE);OUT<-resolve(A$output_dir)
if(normalizePath(OUT,mustWork=FALSE)==DATA)stop("Output and input directories must differ")
source(file.path(ROOT,"assets/style/approved_style.R"))
for(p in c("grid","ragg","svglite","systemfonts"))if(!requireNamespace(p,quietly=TRUE))stop("Missing R package: ",p)
library(grid)
if(!is.null(A$font))STYLE$font<-A$font
if(!STYLE$font%in%systemfonts::system_fonts()$family)stop("Configured font unavailable: ",STYLE$font)
read_tsv<-function(f)read.delim(file.path(DATA,f),check.names=FALSE,stringsAsFactors=FALSE,colClasses="character")
need<-function(d,n){if(!all(n%in%names(d))||!nrow(d)||anyNA(d[,n,drop=FALSE])||any(!nzchar(trimws(as.matrix(d[,n,drop=FALSE])))))stop("Missing required data: ",paste(n,collapse=","))}
num<-function(x){y<-suppressWarnings(as.numeric(x));if(any(!is.finite(y)))stop("Nonfinite numeric data");y}
cells<-read_tsv("cells.tsv");cl<-read_tsv("clusters.tsv")
need(cells,c("cell_id","cluster_id","umap_1","umap_2"));need(cl,c("cluster_id","label","color"))
cells$umap_1<-num(cells$umap_1);cells$umap_2<-num(cells$umap_2)
if(anyDuplicated(cells$cell_id)||anyDuplicated(cl$cluster_id)||!setequal(cells$cluster_id,cl$cluster_id))stop("Cell/cluster identity mismatch")
if(any(!grepl("^#[0-9A-Fa-f]{6}$",cl$color)))stop("Invalid cluster colors")
if(!"group_id"%in%names(cl))cl$group_id<-"all"
if(!"group_label"%in%names(cl))cl$group_label<-cl$group_id
groups<-unique(cl$group_id)
if(anyDuplicated(rle(cl$group_id)$values))stop("Rows belonging to a display group must be contiguous")
if(any(vapply(split(cl$group_label,cl$group_id),function(x)length(unique(x))!=1,logical(1))))stop("Each group needs one label")
cl$n_cells<-as.integer(table(factor(cells$cluster_id,levels=cl$cluster_id)))
mini_mode<-if(is.null(A$mini_mode))"cluster" else A$mini_mode
placement<-if(is.null(A$mini_placement))"strip" else A$mini_placement
stopifnot(mini_mode%in%c("cluster","group","custom","none"),placement%in%c("strip","sidebar"))
if(placement=="sidebar"&&!mini_mode%in%c("group","custom"))stop("Sidebar requires group or custom highlights")
# Membership is independent of display groups. Anchors position sidebar panels only.
hm<-data.frame(highlight_id=character(),label=character(),cluster_id=character(),anchor_group=character())
if(mini_mode=="custom") {
  if(is.null(A$highlights_file))stop("Custom highlights need highlights_file")
  hm<-read_tsv(A$highlights_file);need(hm,c("highlight_id","label","cluster_id"))
  if(anyDuplicated(paste(hm$highlight_id,hm$cluster_id,sep="\r")))stop("Duplicate member within a highlight")
  if(any(!hm$cluster_id%in%cl$cluster_id))stop("Unknown cluster in highlight membership")
  if(!"anchor_group"%in%names(hm))hm$anchor_group<-""
} else if(mini_mode=="cluster") {
  hm<-data.frame(highlight_id=cl$cluster_id,label=cl$cluster_id,cluster_id=cl$cluster_id,anchor_group=cl$group_id)
} else if(mini_mode=="group") {
  hm<-data.frame(highlight_id=cl$group_id,label=cl$group_label,cluster_id=cl$cluster_id,anchor_group=cl$group_id)
}
highlight_ids<-unique(hm$highlight_id)
for(id in highlight_ids) {
  d<-hm[hm$highlight_id==id,]
  if(length(unique(d$label))!=1)stop("Inconsistent label within highlight: ",id)
  if(length(unique(d$anchor_group))!=1||anyNA(d$anchor_group))stop("Inconsistent anchor within highlight: ",id)
}
highlight_labels<-vapply(highlight_ids,function(id)hm$label[match(id,hm$highlight_id)],character(1))
anchors<-vapply(highlight_ids,function(id)hm$anchor_group[match(id,hm$highlight_id)],character(1))
if(placement=="sidebar"&&(any(!anchors%in%groups)||anyDuplicated(anchors)))stop("Each sidebar highlight needs a distinct valid anchor_group")
highlights<-setNames(lapply(highlight_ids,function(id)cl$cluster_id[cl$cluster_id%in%hm$cluster_id[hm$highlight_id==id]]),highlight_ids)
# These exact selected row indices drive both raster rendering and audit tables.
highlight_indices<-lapply(highlights,function(ids)which(cells$cluster_id%in%ids))
layouts<-if(is.null(A$layouts))"portrait" else A$layouts
if(!length(layouts)||!all(layouts%in%c("portrait","wide")))stop("Unknown layout")
if(!length(A$modules))stop("Choose at least one evidence module")
sc<-read.csv(file.path(ROOT,"assets/style/continuous_scales.csv"),stringsAsFactors=FALSE)
width<-function(s,pt)systemfonts::string_width(s,family=STYLE$font,size=pt,res=72)*25.4/72
modules<-lapply(A$modules,function(m){
  if(!m$type%in%c("bubble","heatmap"))stop("Module must be bubble or heatmap")
  if(is.null(m$labels))m$labels<-m$features
  if(!length(m$features)||anyDuplicated(m$features)||length(m$features)!=length(m$labels))stop("Invalid feature labels")
  d<-read_tsv(m$file);need(d,c("cluster_id","feature",if(m$type=="bubble")c("color_value","size_value") else "value"))
  keys<-paste(d$cluster_id,d$feature,sep="\r");expected<-as.vector(outer(cl$cluster_id,m$features,paste,sep="\r"))
  if(anyDuplicated(keys)||!setequal(keys,expected))stop("Module keys mismatch: ",m$title)
  for(k in intersect(c("value","color_value","size_value"),names(d)))d[[k]]<-num(d[[k]])
  if(length(m$limits)!=2||any(!is.finite(m$limits))||diff(m$limits)<=0)stop("Invalid color limits")
  if(is.null(m$ticks))m$ticks<-format(c(m$limits[1],mean(m$limits),m$limits[2]))
  if(length(m$ticks)!=3)stop("Supply three legend labels")
  if(is.null(m$semantics)||!nzchar(m$semantics))stop("State module numeric meaning")
  if(m$type=="bubble") {
    if(is.null(m$size_max)||!is.finite(m$size_max)||m$size_max<=0||any(d$size_value<0|d$size_value>m$size_max))stop("Invalid bubble sizes")
    if(is.null(m$size_breaks))m$size_breaks<-m$size_max*c(.25,.5,1)
    if(length(m$size_breaks)!=3||any(m$size_breaks<0|m$size_breaks>m$size_max))stop("Invalid bubble legend")
  }
  ramp<-if(is.null(m$ramp))if(m$type=="bubble")"or" else "pathway" else m$ramp
  m$colors<-sc$color[sc$scale==ramp];if(!length(m$colors))stop("Unknown color ramp")
  m$data<-d
  m$width<-max(length(m$features)*5.6+3,width(m$title,8)+3,width(m$legend,6.6)+3,26)
  if(m$type=="bubble")m$width<-max(m$width,width(m$size_legend,6.6)+3)
  m
})
if(any(vapply(modules,function(m)is.null(m$title)||is.null(m$legend),logical(1))))stop("Missing module title/legend")
header<-max(22,max(unlist(lapply(modules,function(m)width(m$labels,7.5))))/sqrt(2)+10)
label_width<-max(17,max(width(cl$label,8))+8)
start<-label_width+4
for(j in seq_along(modules)) {
  m<-modules[[j]];m$center<-start+m$width/2
  m$x<-m$center+(seq_along(m$features)-(length(m$features)+1)/2)*5.6
  modules[[j]]<-m;start<-start+m$width+8
}
table_width<-start-6
last<-modules[[length(modules)]]
table_width<-max(table_width,max(last$x+width(last$labels,7.5)/sqrt(2))+3)
bottom<-if(any(vapply(modules,function(m)m$type=="bubble",logical(1))))43 else 20
rh<-STYLE$row_height_mm;gap<-if(length(groups)>1)7 else 0
gh<-vapply(groups,function(g)sum(cl$group_id==g)*rh,numeric(1))
if(placement=="sidebar")gh[groups%in%anchors]<-pmax(gh[groups%in%anchors],34)
rows_height<-sum(gh)+gap*(length(groups)-1)
top<-bottom+rows_height;row_y<-numeric(nrow(cl));group_y<-numeric(length(groups));cursor<-top
for(k in seq_along(groups)) {
  ix<-which(cl$group_id==groups[k]);group_y[k]<-cursor-gh[k]/2
  row_y[ix]<-group_y[k]+((length(ix)+1)/2-seq_along(ix))*rh
  cursor<-cursor-gh[k]-gap
}
table_height<-top+header
side<-if(placement=="sidebar")max(41,max(width(highlight_labels,7.5))+8) else 0
highlight_y<-group_y[match(anchors,groups)]
mini_rows<-ceiling(length(highlights)/6)
strip_h<-if(placement=="strip")mini_rows*40 else 0
main_mm<-if(is.null(A$main_mm))110 else A$main_mm
if(length(main_mm)!=1||!is.finite(main_mm)||main_mm<=0)stop("Invalid main UMAP size")
portrait_dims<-c(max(190,table_width+side+20,main_mm+20),table_height+strip_h+main_mm+29)
wide_dims<-c(table_width+side+main_mm+27,max(main_mm+15,table_height)+strip_h+22)
span<-max(diff(range(cells$umap_1)),diff(range(cells$umap_2)))*STYLE$coordinate_padding
if(span<=0)stop("Degenerate UMAP")
if("--validate-only"%in%args){cat("Validated",nrow(cells),"cells;",nrow(cl),"clusters;",length(groups),"display groups;",length(modules),"modules.\n");quit(save="no")}
dir.create(OUT,recursive=TRUE,showWarnings=FALSE)
colors<-setNames(cl$color,cl$cluster_id)
xmin<-mean(range(cells$umap_1))-span/2;ymin<-mean(range(cells$umap_2))-span/2
ux<-function(x)(x-xmin)/span;uy<-function(y)(y-ymin)/span
mm<-function(x)unit(x,"mm")
txt<-function(s,x,y,size=7.6,just="centre",rot=0,col="#191919")grid.text(s,mm(x),mm(y),just=just,rot=rot,gp=gpar(fontfamily=STYLE$font,fontsize=size,col=col))
rect<-function(x,y,w,h,color)grid.rect(mm(x),mm(y),mm(w),mm(h),gp=gpar(fill=color,col=NA))
circle<-function(x,y,d,color)grid.circle(mm(x),mm(y),r=mm(d/2),gp=gpar(fill=color,col=NA))
vp<-function(x,y,w,h)pushViewport(viewport(x=mm(x),y=mm(y),width=mm(w),height=mm(h),just=c("left","bottom"),clip="off"))
mapcol<-function(v,m){t<-pmax(0,pmin(1,(v-m$limits[1])/diff(m$limits)));m$colors[round(t*(length(m$colors)-1))+1]}
factor<-if(is.null(A$particle_factor))max(.35,min(1.5,sqrt(STYLE$reference_n_cells/nrow(cells)))) else A$particle_factor
stopifnot(is.finite(factor),factor>0)
RNGkind("Mersenne-Twister","Inversion","Rejection");set.seed(STYLE$draw_seed);order_draw<-sample(seq_len(nrow(cells)))
cloud<-function(selected=NULL) {
  cap<-ragg::agg_capture(width=if(is.null(selected))110 else 30,height=if(is.null(selected))110 else 30,units="mm",res=600,background="white");grid.newpage()
  if(!is.null(selected))grid.points(ux(cells$umap_1),uy(cells$umap_2),default.units="npc",pch=16,size=mm(STYLE$mini_point_mm*factor),gp=gpar(col=adjustcolor(STYLE$grey,alpha.f=STYLE$grey_alpha)))
  ix<-if(is.null(selected))order_draw else order_draw[order_draw%in%selected]
  d<-cells[ix,];grid.points(ux(d$umap_1),uy(d$umap_2),default.units="npc",pch=16,size=mm(if(is.null(selected))STYLE$point_mm*factor else STYLE$mini_point_mm*factor),
                            gp=gpar(col=adjustcolor(colors[d$cluster_id],alpha.f=if(is.null(selected))STYLE$point_alpha else STYLE$mini_alpha)))
  img<-cap();dev.off();img
}
labels<-do.call(rbind,lapply(cl$cluster_id,function(id){d<-cells[cells$cluster_id==id,];i<-which.min((d$umap_1-median(d$umap_1))^2+(d$umap_2-median(d$umap_2))^2);data.frame(cluster_id=id,x=d$umap_1[i],y=d$umap_2[i])}))
if(file.exists(file.path(DATA,"label_positions.tsv"))) {
  z<-read_tsv("label_positions.tsv");need(z,c("cluster_id","x","y"));if(anyDuplicated(z$cluster_id)||!setequal(z$cluster_id,cl$cluster_id))stop("Label IDs mismatch")
  z$x<-num(z$x);z$y<-num(z$y);labels<-z[match(cl$cluster_id,z$cluster_id),]
}
hero<-cloud();small<-lapply(highlight_indices,cloud)
main<-function(x,y,size=110) {
  vp(x,y,size,size);grid.raster(hero,width=unit(1,"npc"),height=unit(1,"npc"),interpolate=FALSE)
  for(i in seq_len(nrow(labels))) {
    xx<-ux(labels$x[i])*size;yy<-uy(labels$y[i])*size
    for(a in seq(0,2*pi,length.out=13)[-13])txt(labels$cluster_id[i],xx+.19*cos(a),yy+.19*sin(a),8.6,col="white")
    txt(labels$cluster_id[i],xx,yy,8.6)
  }
  grid.lines(mm(c(3,3,18)),mm(c(18,4,4)),gp=gpar(lwd=.55))
  txt("UMAP-1",11,1,7);txt("UMAP-2",.5,11,7,rot=90);txt(paste0("n = ",format(nrow(cells),big.mark=",")),23,6,7,just="left");popViewport()
}
matrix_panel<-function(x,y) {
  vp(x,y,table_width,table_height)
  for(i in seq_len(nrow(cl))){txt(cl$label[i],label_width-5,row_y[i],8,just="right");circle(label_width-.8,row_y[i],2.4,colors[cl$cluster_id[i]])}
  for(m in modules) {
    txt(m$title,m$center,top+header-4,8)
    for(j in seq_along(m$features))txt(m$labels[j],m$x[j]-1.2,top+3.5,7.5,just="left",rot=45)
    d<-m$data
    for(i in seq_len(nrow(d))) {
      r<-match(d$cluster_id[i],cl$cluster_id);c<-match(d$feature[i],m$features)
      if(m$type=="bubble")circle(m$x[c],row_y[r],4*sqrt(d$size_value[i]/m$size_max),mapcol(d$color_value[i],m))
      else rect(m$x[c],row_y[r],5.615,rh+.015,mapcol(d$value[i],m))
    }
    yy<-if(bottom==43)25 else 9;w<-18;step<-w/length(m$colors)
    for(i in seq_along(m$colors))rect(m$center-w/2+(i-.5)*step,yy,step+.004,2.5,m$colors[i])
    txt(m$legend,m$center,yy+4,6.6)
    for(i in 1:3)txt(m$ticks[i],m$center-9+(i-1)*9,yy-3.6,6.2)
    if(m$type=="bubble") {
      txt(m$size_legend,m$center,12,6.6)
      for(i in 1:3){xx<-m$center-8.4+(i-1)*8.4;circle(xx,7,4*sqrt(m$size_breaks[i]/m$size_max),"#C8C8C8");txt(format(m$size_breaks[i],trim=TRUE),xx,2.5,6.2)}
    }
  }
  popViewport()
}
strip<-function(x,y,w) {
  for(i in seq_along(small)) {
    row<-(i-1)%/%6;ids<-seq.int(row*6+1,min((row+1)*6,length(small)))
    centers<-if(length(ids)==1)w/2 else seq(18,w-18,length.out=length(ids));cx<-x+centers[match(i,ids)];cy<-y+strip_h-21-row*40
    grid.raster(small[[i]],mm(cx),mm(cy),width=mm(29),height=mm(29),interpolate=FALSE)
    txt(highlight_labels[i],cx,cy+17,7.5)
  }
}
draw<-function(layout) {
  dims<-if(layout=="portrait")portrait_dims else wide_dims;grid.newpage()
  txt(A$title,dims[1]/2,dims[2]-5,9)
  if(!is.null(A$subtitle))txt(A$subtitle,dims[1]/2,dims[2]-11,7,col="#666666")
  if(layout=="portrait") {
    mx<-(dims[1]-table_width-side)/2+side;my<-8
    main((dims[1]-main_mm)/2,table_height+strip_h+15,main_mm)
    if(placement=="strip")strip(10,table_height+10,dims[1]-20)
  } else {
    mx<-main_mm+15+side;my<-strip_h+8
    main(7,strip_h+max(8,(table_height-main_mm)/2+8),main_mm)
    if(placement=="strip")strip(10,4,dims[1]-20)
  }
  matrix_panel(mx,my)
  if(placement=="sidebar")for(i in seq_along(small)) {
    grid.raster(small[[i]],mm(mx-side/2),mm(my+highlight_y[i]),width=mm(27),height=mm(27),interpolate=FALSE)
    txt(highlight_labels[i],mx-side/2,my+highlight_y[i]+16,7.5)
  }
  if(length(groups)>1&&placement!="sidebar")for(i in seq_along(groups))txt(cl$group_label[match(groups[i],cl$group_id)],mx-5,my+group_y[i],7,just="right",rot=90)
}
for(layout in layouts) {
  dims<-if(layout=="portrait")portrait_dims else wide_dims;stem<-file.path(OUT,paste0(A$prefix,"_",layout))
  ragg::agg_png(paste0(stem,".png"),width=dims[1],height=dims[2],units="mm",res=450);draw(layout);dev.off()
  if(Sys.info()[["sysname"]]=="Darwin")quartz(type="pdf",file=paste0(stem,".pdf"),width=dims[1]/25.4,height=dims[2]/25.4,family=STYLE$font)
  else cairo_pdf(paste0(stem,".pdf"),width=dims[1]/25.4,height=dims[2]/25.4,family=STYLE$font)
  draw(layout);dev.off()
  svglite::svglite(paste0(stem,".svg"),width=dims[1]/25.4,height=dims[2]/25.4,bg="white");draw(layout);dev.off()
}
write.table(data.frame(cluster_id=cl$cluster_id,group_id=cl$group_id,y_mm=row_y),file.path(OUT,"row_alignment.tsv"),sep="\t",quote=FALSE,row.names=FALSE)
write.table(labels,file.path(OUT,"label_positions.tsv"),sep="\t",quote=FALSE,row.names=FALSE)
if(length(highlights)) {
  hm$n_cells<-cl$n_cells[match(hm$cluster_id,cl$cluster_id)]
  write.table(hm,file.path(OUT,"highlight_membership.tsv"),sep="\t",quote=FALSE,row.names=FALSE)
  hs<-data.frame(highlight_id=highlight_ids,label=highlight_labels,anchor_group=anchors,
                 n_clusters=lengths(highlights),n_unique_cells=lengths(highlight_indices),
                 anchor_y_mm=if(placement=="sidebar")highlight_y else NA_real_)
  write.table(hs,file.path(OUT,"highlight_summary.tsv"),sep="\t",quote=FALSE,row.names=FALSE)
}
capture.output(dput(list(config=A,style=STYLE,particle_factor=factor,portrait_mm=portrait_dims,wide_mm=wide_dims,group_y=group_y,
  highlight_ids=highlight_ids,highlight_y=highlight_y,highlights=highlights,
  input_md5=tools::md5sum(list.files(DATA,full.names=TRUE)))),file=file.path(OUT,"run_record.R"))
capture.output(sessionInfo(),file=file.path(OUT,"sessionInfo.txt"))
cat("Rendered",paste(layouts,collapse=", "),"for",nrow(cells),"cells and",length(modules),"modules.\n")
