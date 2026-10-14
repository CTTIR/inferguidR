# The scaled chi-square family admits a quantile-ratio solution for one tail.
# For two tails the derivative changes sign once, so each side of its minimum
# has at most one root; bracket on the log ratio rather than a fixed scan.
.ig_sensitivity_chi2_variance_one_sample <- function(p, target_name) {
 if(!identical(target_name,"ratio")) stop("variance sensitivity solves ratio",call.=FALSE)
 df<-as.double(p$n_total)-1;alpha<-p$alpha;target<-p$power
 fail<-function(reason) list(status="numerical_failure",reason=reason)
 if(target<alpha) return(list(status="infeasible",reason="Target power must be at least alpha."))
 if(identical(p$tails,"one")) {
  upper<-.ig_g11_qchisq(alpha,df,FALSE)
  denominator<-.ig_g11_qchisq(1-target,df,TRUE)
  value<-if(target==alpha) 1 else upper/denominator
  lower<-.ig_g11_qchisq(alpha,df,TRUE)/.ig_g11_qchisq(target,df,TRUE)
  if(!is.finite(value)||value<1||!is.finite(lower)||lower<0)
   return(fail("The variance-ratio quantile solution is outside the representable numerical domain."))
  return(list(status="ok",value=value,diagnostics=list(
   method="scaled_chi_square_quantile_ratio",branch_rule="ratio at least one; increasing upper-tail branch",
   roots_found=sort(unique(c(lower,value))),other_branch_roots=if(lower!=value) lower else numeric(),
   root_kind=if(target==alpha) "null_effect_boundary" else "quantile_ratio",
   target_reached=TRUE,iterations=0L)))
 }
 lo<-.ig_g11_qchisq(alpha/2,df,TRUE);hi<-.ig_g11_qchisq(alpha/2,df,FALSE)
 if(!is.finite(lo)||!is.finite(hi)||lo<=0||hi<=lo)
  return(fail("The two-tail critical values could not be resolved for variance sensitivity."))
 # Derivative sign is that of df*log(hi/lo) - (hi-lo)/ratio.
 turning<-(hi-lo)/(df*log(hi/lo))
 if(!is.finite(turning)||turning<=0) return(fail("The power-curve minimum could not be resolved."))
 pow<-function(u) {
  ratio<-exp(u)
  stats::pchisq(lo/ratio,df)+stats::pchisq(hi/ratio,df,lower.tail=FALSE)
 }
 f<-function(u) pow(u)-target
 # The selected root is the increasing one, to the right of the minimum.
 centre<-log(turning);fc<-f(centre)
 if(!is.finite(fc)||fc>1e-12) return(fail("The requested power is below the resolved power-curve minimum."))
 solve_side<-function(sign) {
  outside<-centre+sign;fo<-f(outside);iter<-0L
  limit<-if(sign>0) log(.Machine$double.xmax) else log(.Machine$double.xmin)+log(.Machine$double.eps)
  while(is.finite(fo)&&fo<0&&iter<12L&&outside!=limit) {
   span<-abs(outside-centre)*2
   outside<-if(sign>0) min(limit,centre+span) else max(limit,centre-span)
   fo<-f(outside);iter<-iter+1L
  }
  if(!is.finite(fo)||fo<0) return(NULL)
  bracket<-sort(c(centre,outside))
  root<-stats::uniroot(f,bracket,tol=1e-13,maxiter=1000L)
  list(value=exp(root$root),bracket=exp(bracket),iterations=root$iter+iter)
 }
 right<-solve_side(1);left<-solve_side(-1)
 if(is.null(right)||is.null(left)||!is.finite(right$value)||right$value<=0)
  return(fail("A variance-ratio root is outside the representable search domain."))
 list(status="ok",value=right$value,diagnostics=list(
  method="scaled_chi_square_log_ratio_bracketing",branch_rule="increasing branch to the right of the power minimum",
  minimum_ratio=turning,roots_found=c(left$value,right$value),other_branch_roots=left$value,
  bracket=right$bracket,iterations=right$iterations,root_kind="crossing",target_reached=TRUE))
}
