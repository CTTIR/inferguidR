# Synthetic fixed-factor example: prespecify the factorial model and hypotheses.
# Two or three crossed categorical factors; full interactions, no empty cells.
# Type-III hypotheses use sum-to-zero contrasts and unweighted marginal means.
# All term tests are reported, with Holm adjustment across this declared family.
check_settings <- function(alternative,conf.level) {
 stopifnot(identical(alternative,"two.sided"),is.numeric(conf.level),length(conf.level)==1L,
  is.finite(conf.level),conf.level>0,conf.level<1)
}
analyse <- function(data,factor_columns,independent_units=FALSE,equal_variance_justified=FALSE,
 alternative="two.sided",conf.level=.95,complete_cases=FALSE) {
 check_settings(alternative,conf.level)
 stopifnot(isTRUE(independent_units),isTRUE(equal_variance_justified),is.data.frame(data),
  !anyDuplicated(names(data)),is.character(factor_columns),length(factor_columns)%in%c(2L,3L),
  !anyNA(factor_columns),!anyDuplicated(factor_columns),all(nzchar(factor_columns)),
  !any(factor_columns %in% c("id","y")),all(c("id","y",factor_columns)%in%names(data)),
  nrow(data)>=8L,nrow(data)<=10000L,is.numeric(data$y),is.null(dim(data$y)),
  is.atomic(data$id),is.null(dim(data$id)),!anyNA(data$id),!anyDuplicated(data$id),
  all(nzchar(as.character(data$id))),is.logical(complete_cases),length(complete_cases)==1L,!is.na(complete_cases))
 stopifnot(all(vapply(data[factor_columns],function(x)is.factor(x)&&!is.ordered(x)&&
  nlevels(x)>=2L&&!anyNA(levels(x))&&all(nzchar(levels(x))),TRUE)))
 keep<-stats::complete.cases(data[c("y",factor_columns)])
 if(any(!keep)&&!complete_cases)stop("Missing selected measurements require an explicit complete-case choice.")
 n_total<-nrow(data);data<-data[keep,,drop=FALSE];stopifnot(nrow(data)>0,all(is.finite(data$y)))
 levels_count<-vapply(data[factor_columns],nlevels,0L);stopifnot(prod(levels_count)<=64)
 counts<-do.call(table,c(unname(data[factor_columns]),list(dnn=factor_columns)))
 if(any(counts<2L))stop("Every declared crossed-factor cell needs at least two independent observations; no pooling or level dropping is performed.")
 d<-data.frame(.y=data$y-data$y[1]);internal<-paste0(".factor",seq_along(factor_columns))
 for(j in seq_along(internal))d[[internal[j]]]<-data[[factor_columns[j]]]
 contrast_list<-stats::setNames(lapply(levels_count,stats::contr.sum),internal)
 formula<-stats::as.formula(paste(".y ~",paste(internal,collapse=" * ")),env=baseenv())
 fit<-stats::lm(formula,data=d,contrasts=contrast_list,na.action=stats::na.fail,singular.ok=FALSE)
 X<-stats::model.matrix(fit);beta<-stats::coef(fit)
 cell_id<-rep(1,nrow(d));stride<-1
 for(j in seq_along(internal)) {
  cell_id<-cell_id+(as.integer(d[[internal[j]]])-1)*stride;stride<-stride*levels_count[j]
 }
 # In this saturated cell-means model the error is within-cell variation.
 # Compute it directly to avoid tiny QR residuals inventing variance in exact fits.
 df_error<-stats::df.residual(fit)
 sse<-sum(vapply(split(d$.y,cell_id),function(v)sum((v-mean(v))^2),0))
 if(!is.finite(sse)||sse<=0||df_error<=0) return(list(status="not_applicable",p_value=NULL,
  reason="Positive finite residual variation and degrees of freedom are required.",confidence_interval=NULL,interval_status="not_computed"))
 # A quadratic form in the tested coefficient block avoids subtracting two
 # nearly equal residual sums of squares. Sum contrasts define the hypotheses.
 R<-qr.R(fit$qr);V<-chol2inv(R);pivot<-fit$qr$pivot
 unscaled<-matrix(0,ncol(X),ncol(X));unscaled[pivot,pivot]<-V
 assignment<-attr(X,"assign");labels<-attr(stats::terms(fit),"term.labels")
 result<-lapply(seq_along(labels),function(j){
  index<-which(assignment==j);b<-beta[index];q<-length(index)
  ss<-as.numeric(crossprod(b,solve(unscaled[index,index,drop=FALSE],b)))
  F<-(ss/q)/(sse/df_error)
  data.frame(term=labels[j],df=q,sum_squares=ss,F=F,p_value=stats::pf(F,q,df_error,lower.tail=FALSE),
   partial_eta_squared=if(ss==0) 0 else 1/(1+sse/ss),stringsAsFactors=FALSE)
 })
 result<-do.call(rbind,result)
 result$term<-vapply(strsplit(labels,":",fixed=TRUE),function(tokens)
  paste(factor_columns[match(tokens,internal)],collapse=":"),"")
 if(any(!is.finite(as.matrix(result[-1])))||any(result$sum_squares<0))return(list(status="not_applicable",p_value=NULL,
  reason="The factorial F statistics are not numerically resolved.",confidence_interval=NULL,interval_status="not_computed"))
 result$p_holm<-stats::p.adjust(result$p_value,method="holm")
 list(status="ok",method="Fixed-factor full factorial ANOVA; type-III sum-contrast hypotheses",
  terms=result,p_value=stats::setNames(result$p_value,result$term),
  residual_df=df_error,residual_sum_squares=sse,n_total=n_total,n_used=nrow(data),n_excluded=n_total-nrow(data),
  cell_sizes=counts,cell_means=tapply(data$y,data[factor_columns],mean),factor_levels=lapply(data[factor_columns],levels),
  alpha=1-conf.level,confidence_interval=NULL,interval_status="not_computed",requested_confidence_level=conf.level,
  multiplicity="Raw term p-values are marginal; p_holm adjusts the prespecified family of all main effects and interactions.",
  effect_size_note="Type-III sums of squares need not add to the total in unbalanced designs; partial eta squared is not an additive fraction of total explained variation.",
  interpretation="Main effects compare unweighted marginal means. Interactions test departure from the corresponding lower-order effects. Specify simple contrasts and their multiplicity before analysis; none are selected here.",
  assumptions="Independent units, fixed crossed factors, normal errors with a common variance for exact F inference. Complete-case inference additionally needs justified missingness assumptions. No random, repeated, nested or continuous factors.",
  limitations="Two or three factors, at most64 cells and10000 observations; at least two observations in every declared cell. No effect-size sampling interval is implemented.")
}
example_data<-expand.grid(A=factor(c("a1","a2")),B=factor(c("b1","b2","b3")),replicate=1:3)
example_data$id<-seq_len(nrow(example_data))
example_data$y<-2*as.integer(example_data$A)+as.integer(example_data$B)+
 as.integer(example_data$A)*as.integer(example_data$B)+c(-1,0,1)[example_data$replicate]
example_result<-analyse(example_data,c("A","B"),TRUE,TRUE,alternative,confidence_level)
print(example_result)
