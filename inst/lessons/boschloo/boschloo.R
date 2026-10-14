# Independent-binomial teaching example; all counts are synthetic.
run_boschloo_lesson <- function(output_dir = tempfile('inferguidR-boschloo-')) {
  if (!requireNamespace('inferguidR', quietly = TRUE)) stop('Install inferguidR before running this lesson.')
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  table <- matrix(c(7, 1, 1, 7), 2, byrow = TRUE,
    dimnames = list(c('group_1', 'group_2'), c('event', 'non_event')))
  observed <- inferguidR::boschloo_test(7, 8, 1, 8, variant = 'minlike', alpha = .05)
  stopifnot(observed$status$code == 'ok')
  fisher <- stats::fisher.test(table, alternative = 'two.sided')
  marginal_intervals <- rbind(
    group_1 = stats::binom.test(7, 8)$conf.int,
    group_2 = stats::binom.test(1, 8)$conf.int)
  colnames(marginal_intervals) <- c('lower', 'upper')
  # Future probabilities are design assumptions, not the observed proportions.
  plan <- inferguidR::boschloo_sample_size(.8, .2, power = .8, alpha = .05,
    variant = 'minlike', allocation_ratio = 1, stability_window = 3L, n_max = 40L)
  stopifnot(plan$status$code == 'ok', plan$first$power >= .8, plan$minimality$first_is_minimum)
  verified <- inferguidR::boschloo_power(plan$first$n1, plan$first$n2, .8, .2,
    alpha = .05, variant = 'minlike')
  stopifnot(abs(verified$power - plan$first$power) < 1e-12)
  sensitivity <- do.call(rbind,lapply(c(.7,.8,.9),function(p1) {
    answer <- inferguidR::boschloo_power(plan$first$n1,plan$first$n2,p1,.2,alpha=.05,variant='minlike')
    data.frame(p1=p1,p2=.2,power=answer$power)
  }))
  utils::write.csv(plan$staircase,file.path(output_dir,'planning.csv'),row.names=FALSE)
  utils::write.csv(sensitivity,file.path(output_dir,'sensitivity.csv'),row.names=FALSE)
  utils::write.csv(marginal_intervals,file.path(output_dir,'marginal-intervals.csv'))
  for (locale in c('en','de')) {
    snapshot <- inferguidR::create_report_snapshot(
      inferguidR::analysis_envelope(observed$spec,observed),locale=locale,
      title=if(locale=='en') 'Synthetic independent binomial groups' else 'Synthetische unabhaengige Binomialgruppen')
    for(format in c('json','txt')) inferguidR::export_report(snapshot,
      file.path(output_dir,paste0('observed-',locale,'.',format)),format)
    if(inferguidR::report_pdf_capabilities()$available) inferguidR::export_report(snapshot,
      file.path(output_dir,paste0('observed-',locale,'.pdf')),'pdf')
  }
  grDevices::pdf(file.path(output_dir,'planning.pdf'),width=7,height=5)
  tryCatch({
    graphics::plot(plan$staircase$n_total,plan$staircase$power,pch=19,
      xlab='Total sample size / Gesamtstichprobenumfang',ylab='Power',ylim=c(0,1),
      main='Evaluated designs / Berechnete Entwuerfe')
    graphics::abline(h=.8,lty=2)
    graphics::mtext('Fixed p1=0.8, p2=0.2, alpha=0.05; minlike; allocation 1:1',side=3,cex=.7)
  },finally=grDevices::dev.off())
  summary <- list(table=table,observed_risk_difference=7/8-1/8,
    boschloo_p=observed$p_value,numerical_upper_bound=observed$p_upper,
    fisher_conditional_p=fisher$p.value,
    marginal_intervals=marginal_intervals,
    interval_note='Two unadjusted marginal Clopper-Pearson intervals; neither is an interval for the risk difference.',
    first_design=plan$first,stable_design=plan$stable,sensitivity=sensitivity,
    versions=c(inferguidR=as.character(utils::packageVersion('inferguidR')),R=as.character(getRversion())),
    certification='Development teaching example; full release certification remains pending.')
  capture.output(summary,file=file.path(output_dir,'results.txt'))
  capture.output(utils::sessionInfo(),file=file.path(output_dir,'session-info.txt'))
  invisible(summary)
}
