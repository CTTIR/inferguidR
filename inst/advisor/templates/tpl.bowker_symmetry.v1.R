# Two matched occasions, same named categories in the same order on both axes.
# Cells count participants; each participant contributes once to a confirmed pair.
# Independent participants, complete observed pairs, chi-squared asymptotic inference.
# Missing pairs must be handled and reported before constructing this aggregate table.
check_table <- function(counts, alternative, conf.level) {
 stopifnot(is.matrix(counts),is.numeric(counts),nrow(counts)>=3L,nrow(counts)==ncol(counts),
   all(is.finite(counts)),all(counts>=0),all(counts==floor(counts)),sum(counts)>0,sum(counts)<=2^31-1,
   !is.null(rownames(counts)),identical(rownames(counts),colnames(counts)),
   !anyNA(rownames(counts)),all(nzchar(rownames(counts))),!anyDuplicated(rownames(counts)),
   identical(alternative,"two.sided"),is.numeric(conf.level),length(conf.level)==1L,
   is.finite(conf.level),conf.level>0,conf.level<1)
}
analyse <- function(counts, alternative="two.sided", conf.level=.95, n_excluded=0L) {
 check_table(counts,alternative,conf.level);k<-nrow(counts);n<-sum(counts)
 stopifnot(is.numeric(n_excluded),length(n_excluded)==1L,is.finite(n_excluded),
   n_excluded>=0,n_excluded==floor(n_excluded))
 # Symmetry is stronger than equal margins: p[i,j] = p[j,i].
 pair_totals<-(counts+t(counts))[upper.tri(counts)]
 if(any(pair_totals==0)) return(list(status="unavailable",p_value=NULL,
   reason="A category-pair total is zero; this example does not infer structural-zero degrees of freedom."))
 fit<-stats::mcnemar.test(counts,correct=FALSE)
 statistic<-unname(fit$statistic);df<-unname(fit$parameter);p<-fit$p.value
 null<-"Joint distribution symmetric about its diagonal"
 list(method="bowker_symmetry",status="ok",n_pairs=n,n_excluded=n_excluded,categories=rownames(counts),null=null,
   observed=counts,row_proportions=rowSums(counts)/n,column_proportions=colSums(counts)/n,
   statistic=statistic,df=df,p_value=p,confidence_interval=NULL,interval_status="not_computed",
   requested_confidence_level=conf.level,alpha=1-conf.level,reject=p<1-conf.level,
   assumptions="Independent matched participants; identical category meaning/order; chi-squared large-sample approximation. Symmetry and marginal homogeneity are different hypotheses; no automatic substitution.")
}
example_counts<-matrix(c(20,4,6,8,25,3,5,7,18),3,byrow=TRUE,
 dimnames=list(c("A","B","C"),c("A","B","C")))
example_result<-analyse(example_counts,alternative,confidence_level)
print(example_result)
