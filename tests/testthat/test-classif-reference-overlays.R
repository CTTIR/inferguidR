test_that('ROC reference curves satisfy Bayes at ordinary and extreme prevalence', {
 for(pi in c(1e-9,.05,.5,1-1e-9)) {
  r<-evaluate_cutoff(update(classification_examples('scored_synthetic')$spec,prevalence=pi))
  pd<-classification_plot_data(r,'roc',ppv_levels=c(.1,.5,.9));d<-pd$data$iso_ppv
  good<-d$fpr>0
  posterior<-pi*d$tpr[good]/(pi*d$tpr[good]+(1-pi)*d$fpr[good])
  expect_equal(posterior,d$ppv_level[good],tolerance=1e-12)
  expect_true(all(is.finite(d$tpr)&is.finite(d$fpr)))
  expect_true(all(d$tpr>=0&d$tpr<=1&d$fpr>=0&d$fpr<=1))
  for(level in c(.1,.5,.9)) {
   z<-d[d$ppv_level==level,];expect_gt(nrow(z),1)
   expect_equal(max(c(z$tpr,z$fpr)),1,tolerance=1e-12)
  }
 }
})
test_that('F1 references meet the unit square and retain exact harmonic means', {
 r<-evaluate_cutoff(classification_examples('scored_synthetic')$spec)
 d<-classification_plot_data(r,'pr',f1_levels=c(.01,.2,.8,.99))$data$iso_f1
 expect_equal(2*d$recall*d$precision/(d$recall+d$precision),d$f1_level,tolerance=1e-12)
 for(level in c(.01,.2,.8,.99)) {
  z<-d[d$f1_level==level,];expect_equal(z$precision[1],1,tolerance=1e-12)
  expect_equal(tail(z$recall,1),1,tolerance=1e-12)
 }
})
test_that('rendered mathematical references do not join attained PR points', {
 skip_if_not_installed('ggplot2')
 r<-evaluate_cutoff(update(classification_examples('scored_synthetic')$spec,prevalence=.05))
 rr<-data.frame(id=c('variable','empty'),tp=c(3,0),fn=c(2,0),fp=c(1,0),tn=c(4,0))
 for(lang in c('en','de')) for(panel in c('roc','pr')) {
  # The zero-total rule is invalid rather than silently plotted at the origin.
  expect_error(plot_classification(r,panel,realised_rules=rr,lang=lang),class='inferguidR_invalid_input')
  rules<-rr[1,,drop=FALSE];p<-plot_classification(r,panel,realised_rules=rules,lang=lang)
  expect_no_warning(ggplot2::ggplot_build(p));pd<-attr(p,'inferguidR_plot_data')
  expect_equal(nrow(pd$data$realised_rules),1)
  paths<-Filter(function(l) inherits(l$geom,'GeomPath'),p$layers)
  if(panel=='pr')expect_true(all(vapply(paths,function(l)'f1_level'%in%names(l$data),TRUE)))
  triangles<-Filter(function(l)identical(l$aes_params$shape,17),p$layers)
  expect_length(triangles,1);expect_equal(triangles[[1]]$data$marker,'realised_rule')
 }
 p<-plot_classification(r,'pr',basis='scenario',realised_rules=rr[1,,drop=FALSE])
 expect_null(attr(p,'inferguidR_plot_data')$data$realised_rules)
 expect_match(p$labels$caption,'not shown')
 for(panel in c('mosaic','metric_network','metric_vs_cutoff'))
  expect_error(plot_classification(r,panel,realised_rules=rr[1,,drop=FALSE]),class='inferguidR_invalid_input')
 p<-plot_classification(r,'metric_pair',realised_rules=rr[1,,drop=FALSE])
 expect_no_warning(ggplot2::ggplot_build(p))
 expect_equal(attr(p,'inferguidR_plot_data')$data$realised_rules$x,.2)
 expect_equal(attr(p,'inferguidR_plot_data')$data$realised_rules$y,.6)
 p<-plot_classification(r,'metric_pair',metric_x='ppv',basis='scenario',realised_rules=rr[1,,drop=FALSE])
 expect_null(attr(p,'inferguidR_plot_data')$data$realised_rules)
 expect_match(p$labels$caption,'hidden')
})
