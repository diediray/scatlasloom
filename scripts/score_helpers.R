# Expression and pathway summaries; biological choices stay in project scripts.
symbol_logcounts <- function(counts, symbols=rownames(counts), requested=unique(symbols)) {
  stopifnot(inherits(counts,"sparseMatrix"),length(symbols)==nrow(counts),!anyNA(symbols),
            all(counts@x>=0),all(counts@x==floor(counts@x)),!anyDuplicated(colnames(counts)))
  totals<-Matrix::colSums(counts);stopifnot(all(totals>0))
  present<-intersect(requested,unique(symbols));ix<-which(symbols %in% present)
  collapse<-Matrix::sparseMatrix(i=match(symbols[ix],present),j=seq_along(ix),x=1,
                               dims=c(length(present),length(ix)))
  raw<-collapse%*%counts[ix,,drop=FALSE];rownames(raw)<-present;colnames(raw)<-colnames(counts)
  norm<-raw%*%Matrix::Diagonal(x=10000/totals);norm@x<-log1p(norm@x)
  list(raw=raw,logcounts=norm,missing=setdiff(requested,present))
}
cluster_summary <- function(x,cluster,levels=unique(cluster)) {
  stopifnot(length(cluster)==ncol(x))
  out<-matrix(vapply(levels,function(k)Matrix::rowMeans(x[,cluster==k,drop=FALSE]),numeric(nrow(x))),nrow=nrow(x),ncol=length(levels))
  dimnames(out)<-list(rownames(x),levels);out
}
column_z <- function(x) {
  y<-x;for(j in seq_len(ncol(x))) {s<-sd(x[,j]);y[,j]<-if(is.finite(s)&&s>0)(x[,j]-mean(x[,j]))/s else 0};y
}
column_relative <- function(x) {
  y<-x;for(j in seq_len(ncol(x))) {d<-diff(range(x[,j]));y[,j]<-if(d>0)(x[,j]-min(x[,j]))/d else 0};y
}
program_means <- function(gene_means,programs) {
  out<-vapply(programs,function(genes) {
    hit<-intersect(genes,rownames(gene_means));if(!length(hit))stop("Program has no measured genes")
    colMeans(gene_means[hit,,drop=FALSE])
  },numeric(ncol(gene_means)))
  out<-matrix(out,nrow=ncol(gene_means),ncol=length(programs),dimnames=list(colnames(gene_means),names(programs)));out
}
progeny_mlm <- function(logcounts,network,top=100L,block_size=256L) {
  required<-c("gene","pathway","weight","p.value");stopifnot(all(required%in%names(network)))
  net<-do.call(rbind,lapply(split(network,network$pathway),function(d)d[head(order(d$p.value,d$gene),top),]))
  net$measured<-net$gene%in%rownames(logcounts)
  genes<-sort(unique(net$gene[net$measured]));paths<-sort(unique(net$pathway))
  W<-matrix(0,length(genes),length(paths),dimnames=list(genes,paths));d<-net[net$measured,]
  stopifnot(!anyDuplicated(paste(d$gene,d$pathway)))
  W[cbind(match(d$gene,genes),match(d$pathway,paths))]<-d$weight
  X<-cbind(intercept=1,W);q<-qr(X);stopifnot(q$rank==ncol(X),nrow(X)>ncol(X))
  iv<-chol2inv(chol(crossprod(X)));df<-nrow(X)-ncol(X)
  result<-matrix(NA_real_,length(paths),ncol(logcounts),dimnames=list(paths,colnames(logcounts)))
  for(start in seq.int(1L,ncol(logcounts),by=block_size)) {
    ix<-start:min(start+block_size-1L,ncol(logcounts));Y<-as.matrix(logcounts[genes,ix,drop=FALSE])
    B<-qr.coef(q,Y);E<-Y-X%*%B;s2<-colSums(E*E)/df
    result[,ix]<-(B/sqrt(outer(diag(iv),pmax(s2,.Machine$double.eps))))[-1,,drop=FALSE]
    if(start==1L) {
      error<-max(abs(summary(lm(Y[,1]~W))$coefficients[-1,3]-result[,ix[1]]))
      stopifnot(error<1e-8)
    }
  }
  stopifnot(all(is.finite(result)))
  list(activity=result,coverage=net,lm_check_max_error=error)
}

# Exact paired sign-flip test for mean patient-level percentage-point changes.
paired_composition <- function(metadata,cluster_levels,min_cells=1L) {
  stopifnot(all(c("patient_id","condition","cluster_id")%in%names(metadata)))
  stopifnot(all(metadata$condition%in%c("Normal","Tumor")))
  counts<-table(factor(metadata$patient_id),factor(metadata$condition,levels=c("Normal","Tumor")),
                factor(metadata$cluster_id,levels=cluster_levels))
  denom<-apply(counts,c(1,2),sum);paired<-rownames(denom)[apply(denom>=min_cells,1,all)]
  if(length(paired)<2)stop("Insufficient paired patients")
  if(length(paired)>20)stop("Exact enumeration would exceed 2^20; choose a project-appropriate alternative")
  prop<-sweep(counts,c(1,2),pmax(denom,1),"/")*100
  delta<-prop[paired,"Tumor",,drop=FALSE]-prop[paired,"Normal",,drop=FALSE]
  dim(delta)<-c(length(paired),length(cluster_levels));dimnames(delta)<-list(paired,cluster_levels)
  signs<-as.matrix(expand.grid(rep(list(c(-1,1)),length(paired))))
  p<-vapply(seq_along(cluster_levels),function(j)mean(abs(as.vector(signs%*%delta[,j])/length(paired))>=abs(mean(delta[,j]))-1e-12),numeric(1))
  list(summary=data.frame(cluster_id=cluster_levels,delta_pp=colMeans(delta),p=p,p_adj=p.adjust(p,"BH"),n_pairs=length(paired)),
       delta=delta,counts=counts,denominator=denom,patients=paired,min_cells=min_cells)
}
