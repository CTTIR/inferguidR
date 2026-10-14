# Invert the increasing noncentrality parameter, then convert lambda=N*w^2.
.ig_sensitivity_chi2_goodness_of_fit <- function(p, target_name) {
 if(!identical(target_name,"w")) stop("goodness-of-fit sensitivity solves w",call.=FALSE)
 if(p$power<p$alpha) return(list(status="infeasible",reason="Target power must be at least alpha."))
 if(p$power==p$alpha) return(list(status="ok",value=0,diagnostics=list(
  method="noncentral_chi_square_boundary",root_kind="null_effect_boundary",target_reached=TRUE,iterations=0L)))
 at<-function(lambda) {
  q<-p;q$w<-sqrt(lambda/p$n_total)
  .ig_power_chi2_goodness_of_fit_contingency(q)$power-p$power
 }
 low<-0;fl<-at(low)
 if(!is.finite(fl)||fl>=0||abs(fl)<=.ig_solver_settings$power_noise)
  return(list(status="numerical_failure",reason="Target power cannot be separated from the null-effect boundary at the forward model's precision."))
 high<-1;fh<-at(high);iterations<-0L
 while(is.finite(fh)&&fh<0&&iterations<1023L) {
  low<-high;fl<-fh
  if(high>.Machine$double.xmax/2) break
  high<-high*2;fh<-at(high);iterations<-iterations+1L
 }
 if(!is.finite(fh)||fh<0) return(list(status="numerical_failure",
  reason="The noncentrality needed for the target is outside the forward model's resolved computational domain.",
  diagnostics=list(noncentrality_bracket=c(low,high),iterations=iterations)))
 root<-.ig_refine_root(at,low,high,fl,fh)
 if(!root$kind%in%c("crossing","steep_crossing")) return(list(status="numerical_failure",
  reason="The continuous noncentrality root could not be resolved.",diagnostics=list(root_kind=root$kind)))
 value<-sqrt(root$x/p$n_total)
 if(!is.finite(value)||value<0) return(list(status="numerical_failure",reason="The solved effect is not representable."))
 list(status="ok",value=value,diagnostics=list(method="noncentral_chi_square_monotone_bracketing",
  noncentrality=root$x,noncentrality_bracket=root$bracket,iterations=iterations+root$iterations,
  root_kind=root$kind,target_reached=TRUE))
}
