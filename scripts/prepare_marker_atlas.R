#!/usr/bin/env Rscript
# Configurable marker evidence preparation. No clustering, relabeling or pathway analysis.
suppressPackageStartupMessages(library(Matrix))
args<-commandArgs(TRUE);if(!length(args))stop('Usage: prepare_marker_atlas.R config.R')
config<-normalizePath(args[1]);e<-new.env(parent=baseenv());sys.source(config,e);J<-e$MARKER_JOB
if(is.null(J))stop('Config must define MARKER_JOB')
resolve<-function(p)if(grepl('^(/|[A-Za-z]:[/\\\\])',p))p else file.path(dirname(config),p)
raw_path<-normalizePath(resolve(J$object_path));marker_path<-normalizePath(resolve(J$markers_path))
OUT<-resolve(J$output_dir);for(d in c('input','analysis','provenance'))dir.create(file.path(OUT,d),recursive=TRUE,showWarnings=FALSE)
write_tsv<-function(x,n,folder='analysis')write.table(x,file.path(OUT,folder,n),sep='\t',quote=FALSE,row.names=FALSE)
before<-tools::md5sum(c(raw_path,marker_path));o<-readRDS(raw_path)
if(inherits(o,'Seurat')) {
 if(!requireNamespace('SeuratObject',quietly=TRUE))stop('SeuratObject is needed to read this input')
 assay<-if(is.null(J$assay))'RNA' else J$assay
 counts<-SeuratObject::LayerData(o[[assay]],layer='counts');meta<-o[[]]
 reduction<-if(is.null(J$reduction))'umap' else J$reduction
 coords<-SeuratObject::Embeddings(o,reduction)
} else if(is.list(o)&&all(c('counts','metadata','umap')%in%names(o))) {
 counts<-o$counts;meta<-o$metadata;coords<-as.matrix(o$umap)
} else stop('Use a Seurat RDS, or an adapter RDS list containing counts, metadata and umap')
stopifnot(inherits(counts,'sparseMatrix'),!is.null(rownames(counts)),!is.null(colnames(counts)),
          !anyDuplicated(rownames(counts)),!anyDuplicated(colnames(counts)),all(counts@x>=0),all(counts@x==floor(counts@x)),
          !anyDuplicated(rownames(meta)),setequal(rownames(meta),colnames(counts)),
          !anyDuplicated(rownames(coords)),setequal(rownames(coords),colnames(counts)),ncol(coords)==2)
meta<-meta[colnames(counts),,drop=FALSE];coords<-coords[colnames(counts),,drop=FALSE]
if(!J$cluster_field%in%names(meta))stop('Configured cluster_field is missing')
cluster<-as.character(meta[[J$cluster_field]]);cls<-J$cluster_order
if(is.null(cls))cls<-unique(cluster)
stopifnot(length(cls)>=2,!anyNA(cluster),all(nzchar(cluster)),!anyDuplicated(cls),setequal(cluster,cls))
if(is.null(J$colors)||!setequal(names(J$colors),cls)||any(!grepl('^#[0-9A-Fa-f]{6}$',J$colors)))stop('Provide named colors for all clusters')
coordinate_error<-0
if(!is.null(J$coordinates_path)) {
 cp<-normalizePath(resolve(J$coordinates_path));z<-read.delim(cp,check.names=FALSE,stringsAsFactors=FALSE)
 map<-J$coordinate_columns;stopifnot(all(c('cell','x','y')%in%names(map)),all(map%in%names(z)),!anyDuplicated(z[[map[['cell']]]]))
 ix<-match(colnames(counts),z[[map[['cell']]]]);stopifnot(!anyNA(ix))
 frozen<-as.matrix(z[ix,map[c('x','y')],drop=FALSE]);coordinate_error<-max(abs(coords-frozen))
 if(coordinate_error>1e-10)stop('External coordinates disagree with object embedding; inspect the intended frozen source')
 coords<-frozen;before<-c(before,tools::md5sum(cp))
}
stopifnot(all(is.finite(coords)))
total<-Matrix::colSums(counts);stopifnot(all(total>0))
linear<-counts%*%Diagonal(x=10000/total);dimnames(linear)<-dimnames(counts)
logrna<-linear;logrna@x<-log1p(logrna@x)
genes<-rownames(counts);mean_linear<-vapply(cls,function(k)Matrix::rowMeans(linear[,cluster==k,drop=FALSE]),numeric(nrow(linear)))
mean_log<-vapply(cls,function(k)Matrix::rowMeans(logrna[,cluster==k,drop=FALSE]),numeric(nrow(logrna)))
pct<-vapply(cls,function(k)100*Matrix::rowMeans(counts[,cluster==k,drop=FALSE]>0),numeric(nrow(counts)))
dimnames(mean_linear)<-dimnames(mean_log)<-dimnames(pct)<-list(genes,cls)
con<-if(grepl('\\.gz$',marker_path))gzfile(marker_path,open='rt') else file(marker_path,open='rt')
m<-read.delim(con,check.names=FALSE,stringsAsFactors=FALSE,colClasses='character');close(con)
columns<-J$marker_columns
if(!all(c('group','feature','padj')%in%names(columns))||!all(columns%in%names(m)))stop('Map marker group, feature and adjusted-P columns; optionally AUC')
for(k in names(columns))m[[k]]<-m[[columns[[k]]]]
m$group<-as.character(m$group);m$feature<-as.character(m$feature)
m$padj<-suppressWarnings(as.numeric(m$padj))
stopifnot(all(m$group%in%cls),all(cls%in%m$group),all(m$feature%in%genes),!anyDuplicated(paste(m$group,m$feature,sep='\r')),
          is.numeric(m$padj),all(is.finite(m$padj)),all(m$padj>=0&m$padj<=1))
has_auc<-'auc'%in%names(columns)
if(has_auc)m$auc<-suppressWarnings(as.numeric(m$auc))
if(has_auc)stopifnot(is.numeric(m$auc),all(is.finite(m$auc)),all(m$auc>=0&m$auc<=1))
F<-J$filters;needed<-c('max_adjusted_p','min_pct_in','min_delta_pct','min_log2fc','max_per_cluster')
if(!all(needed%in%names(F)))stop('Specify display filters in the project config')
if(!is.null(F$min_auc)&&!has_auc)stop('AUC filter configured, but no AUC column supplied')
stopifnot(F$max_per_cluster>=1,F$max_per_cluster==as.integer(F$max_per_cluster),F$min_pct_in>=0,F$min_pct_in<=100)
sizes<-as.integer(table(factor(cluster,levels=cls)));total_linear<-as.vector(mean_linear%*%sizes)
total_detected<-as.vector(pct%*%sizes)
g<-match(m$feature,genes);k<-match(m$group,cls);inside<-cbind(g,k)
m$avg_log2FC<-log2((mean_linear[inside]+1)/((total_linear[g]-mean_linear[inside]*sizes[k])/(sum(sizes)-sizes[k])+1))
if('pct_in'%in%names(m))m$source_pct_in<-m$pct_in
if('pct_out'%in%names(m))m$source_pct_out<-m$pct_out
m$pct_in<-pct[inside];m$pct_out<-(total_detected[g]-pct[inside]*sizes[k])/(sum(sizes)-sizes[k])
m$delta_pct_pp<-m$pct_in-m$pct_out
m$technical_family<-if(is.null(F$exclude_regex)||!nzchar(F$exclude_regex))FALSE else grepl(F$exclude_regex,m$feature)
m$pass_display_filter<-m$padj<F$max_adjusted_p&m$pct_in>=F$min_pct_in&m$delta_pct_pp>=F$min_delta_pct&m$avg_log2FC>=F$min_log2fc&!m$technical_family
if(!is.null(F$min_auc))m$pass_display_filter<-m$pass_display_filter&m$auc>=F$min_auc
m$selection_reason<-ifelse(m$pass_display_filter,'eligible','below project display criteria or excluded family')
ord<-if(has_auc)order(match(m$group,cls),-m$auc,-m$avg_log2FC,-m$delta_pct_pp,m$feature) else order(match(m$group,cls),-m$avg_log2FC,-m$delta_pct_pp,m$feature)
m<-m[ord,]
policy<-if(is.null(J$donor_policy))'report' else J$donor_policy
stopifnot(policy%in%c('none','report','filter'))
patient<-if(is.null(J$patient_field))rep('',nrow(meta)) else as.character(meta[[J$patient_field]])
if(!is.null(J$patient_field))stopifnot(J$patient_field%in%names(meta),length(patient)==nrow(meta),!anyNA(patient),all(nzchar(patient)))
if(policy=='filter'&&is.null(J$patient_field))stop('Donor filtering needs an identified biological donor field')
m$informative_patients<-0L;m$patient_fraction_positive<-NA_real_;m$patient_consistent<-NA
m$pass_final<-m$pass_display_filter
if(policy!='none'&&!is.null(J$patient_field)) {
 D<-J$donor_rules;if(!all(c('min_cells_each_side','min_patients','min_positive_fraction')%in%names(D)))stop('Specify donor_rules')
 support<-do.call(rbind,lapply(cls,function(k){
  gs<-m$feature[m$group==k&m$pass_display_filter];if(!length(gs))return(NULL)
  do.call(rbind,lapply(sort(unique(patient)),function(p){
   a<-which(patient==p&cluster==k);b<-which(patient==p&cluster!=k)
   if(length(a)<D$min_cells_each_side||length(b)<D$min_cells_each_side)return(NULL)
   aa<-Matrix::rowMeans(linear[gs,a,drop=FALSE]);bb<-Matrix::rowMeans(linear[gs,b,drop=FALSE])
   data.frame(cluster_id=k,gene=gs,patient_id=p,n_in=length(a),n_out=length(b),log2FC=log2((aa+1)/(bb+1)))
  }))
 }))
 if(!is.null(support)&&nrow(support)) {
  write_tsv(support,'within_patient_marker_effects.tsv')
  sm<-do.call(rbind,lapply(split(support,paste(support$cluster_id,support$gene,sep=':')),function(d)data.frame(cluster_id=d$cluster_id[1],gene=d$gene[1],informative_patients=nrow(d),positive_patients=sum(d$log2FC>0),median_patient_log2FC=median(d$log2FC))))
  sm$patient_fraction_positive<-sm$positive_patients/sm$informative_patients
  sm$patient_consistent<-sm$patient_fraction_positive>=D$min_positive_fraction&sm$median_patient_log2FC>0
  write_tsv(sm,'marker_patient_support.tsv')
  ix<-match(paste(m$group,m$feature,sep=':'),paste(sm$cluster_id,sm$gene,sep=':'));hit<-!is.na(ix)
  m$informative_patients[hit]<-sm$informative_patients[ix[hit]]
  m$patient_fraction_positive[hit]<-sm$patient_fraction_positive[ix[hit]]
  m$patient_consistent[hit]<-sm$patient_consistent[ix[hit]]
 }
 assessable<-m$informative_patients>=D$min_patients
 if(policy=='filter')m$pass_final<-m$pass_display_filter&(!assessable|(!is.na(m$patient_consistent)&m$patient_consistent))
 m$selection_reason[m$pass_display_filter&!m$pass_final]<-'pooled candidate; below configured donor consistency'
 m$selection_reason[m$pass_final&!assessable]<-'pooled descriptive candidate; too few informative donors'
}
selected<-do.call(rbind,lapply(cls,function(k)head(m[m$group==k&m$pass_final,],F$max_per_cluster)))
write_tsv(m,'all_marker_selection_evidence.tsv')
selection_counts<-data.frame(cluster_id=cls,n_eligible=sapply(cls,function(k)sum(m$group==k&m$pass_display_filter)),n_selected=sapply(cls,function(k)sum(selected$group==k)))
write_tsv(selection_counts,'selection_counts.tsv')
if(any(selection_counts$n_selected==0))stop('Some clusters have no qualifying marker. Inspect saved candidate evidence before revising the display design')
selected$marker_rank<-ave(seq_len(nrow(selected)),selected$group,FUN=seq_along);selected$row_id<-paste(selected$group,selected$feature,sep=':')
write_tsv(selected,'selected_markers.tsv')
gs<-selected$feature;dot_mean_z<-t(scale(t(mean_log[gs,,drop=FALSE])))
cell_expr<-as.matrix(logrna[gs,,drop=FALSE]);cell_z<-(cell_expr-rowMeans(cell_expr))/apply(cell_expr,1,sd)
stopifnot(all(is.finite(dot_mean_z)),all(is.finite(cell_z)))
cells<-data.frame(cell_id=colnames(counts),cluster_id=cluster,umap_1=coords[,1],umap_2=coords[,2],patient_id=patient)
bundle<-list(title=J$title,cells=cells,clusters=cls,colors=J$colors[cls],markers=selected,dot_mean_z=dot_mean_z,dot_pct=pct[gs,,drop=FALSE],cell_z=cell_z,
             metadata=list(expression='log1p(CP10K)',dot='gene-wise z of cluster mean log expression; size=raw UMI-positive %',
                           heatmap='gene-wise z of log expression across all cells, display clipped at +/-2',selection=F,donor_policy=policy,donor_rules=J$donor_rules,
                           marker_table_scope=J$marker_table_scope,adjusted_p_definition=J$adjusted_p_definition,
                           cell_order='cluster, patient_id, cell_id; equal-width cluster blocks contain all cells'))
saveRDS(bundle,file.path(OUT,'input/marker_atlas.rds'));write_tsv(cells,'cells.tsv','input')
stopifnot(identical(before,tools::md5sum(names(before))))
write_tsv(data.frame(path=names(before),md5=unname(before)),'source_files.tsv','provenance')
capture.output(dput(J),file=file.path(OUT,'provenance/effective_job.R'));capture.output(sessionInfo(),file=file.path(OUT,'provenance/sessionInfo.txt'))
write_tsv(data.frame(check=c('n_cells','n_clusters','n_marker_rows','object_vs_frozen_coordinate_error','original_inputs_unchanged'),value=c(nrow(cells),length(cls),nrow(selected),coordinate_error,TRUE)),'checks.tsv','provenance')
cat('Prepared',nrow(selected),'marker rows for',length(cls),'clusters and',nrow(cells),'cells.\n')
