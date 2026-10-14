# Portable power sweeps retain individual specifications and results, including failures.
.ig_report_sweeps <- function(sweeps,envelopes) {
 bad<-function(message).ig_abort('invalid_input',paste('Invalid report sweep:',message))
 if(!is.list(sweeps)||length(sweeps)>20L)bad('provide a list of at most20 sweep entries')
 lapply(sweeps,function(entry){
  if(!is.list(entry)||!setequal(names(entry),c('analysis','sweep'))||anyDuplicated(names(entry)))bad('entry needs analysis and sweep')
  i<-entry$analysis
  if(!is.numeric(i)||length(i)!=1L||!is.finite(i)||i!=floor(i)||i<1||i>length(envelopes))bad('unknown analysis')
  e<-envelopes[[i]];sw<-entry$sweep;free<-attr(sw,'free')
  if(e$analysis_kind!='power'||!identical(e$result$status$code,'ok'))bad('base analysis must be completed power analysis')
  if(!is.data.frame(sw)||nrow(sw)<1L||nrow(sw)>1000L||!is.numeric(sw$value)||any(!is.finite(sw$value)))bad('sweep requires1\u20131000 finite requested input values')
  if(!.ig_is_string(free)||!free%in%names(e$spec$inputs)||free%in%e$spec$target||
    !identical(attr(sw,'base_spec_hash'),e$spec$spec_hash)||!identical(attr(sw,'method_id'),e$spec$method_id)||
    !identical(attr(sw,'mode'),e$spec$mode))bad('sweep does not belong to the frozen base specification')
  rows<-lapply(seq_len(nrow(sw)),function(j){
   spec<-tryCatch(stats::update(e$spec,inputs=stats::setNames(list(sw$value[j]),free)),inferguidR_invalid_input=function(err)err)
   result<-sw$result[[j]]
   if(inherits(spec,'condition')){
    if(!identical(sw$admissible[j],FALSE)||!is.na(sw$spec_hash[j])||!is.na(sw$result_hash[j])||!is.null(sw$spec[[j]])||!is.null(result)||!identical(sw$status[j],'invalid_input')||!identical(sw$reason[j],conditionMessage(spec)))bad('invalid row does not match its requested input')
    return(list(value=sw$value[j],status='invalid_input',reason=conditionMessage(spec),specification=NULL,result=NULL))
   }
   if(!identical(sw$admissible[j],TRUE)||!inherits(sw$spec[[j]],'power_spec')||!inherits(result,'power_result'))bad('missing typed row specification/result')
   rebuilt<-analysis_envelope(sw$spec[[j]],result)
   if(!identical(spec$spec_hash,rebuilt$spec_hash)||!identical(sw$spec_hash[j],spec$spec_hash)||
     !identical(sw$result_hash[j],result$result_hash)||!identical(sw$status[j],result$status$code)||!identical(sw$reason[j],result$status$reason))bad('row hashes or status disagree')
   number<-function(v)if(!identical(result$status$code,'ok')||is.null(v)||length(v)!=1L)NA_real_ else as.double(v)
   expected<-c(stats::setNames(lapply(e$spec$target,function(t)number(result$solved[[t]])),paste0('solved_',e$spec$target)),
    list(achieved_power=number(result$achieved$power),achieved_beta=number(result$achieved$beta),alpha_attained=number(result$achieved$alpha_attained),n_total=number(result$solved$n_total%||%result$quantities$n_total)))
   for(name in names(expected))if(!name%in%names(sw)||!identical(as.double(sw[[name]][j]),expected[[name]]))bad('displayed sweep values disagree with retained result')
   list(value=sw$value[j],status=result$status$code,reason=result$status$reason,specification=.ig_spec_to_list(spec),result=.ig_result_to_list(result))
  })
  list(analysis=as.integer(i),schema_version='inferguidR.report_power_sweep/0.1.0',base_spec_hash=e$spec$spec_hash,free=free,rows=rows)
 })
}

.ig_report_validate_sweeps <- function(x) {
 bad<-function(message).ig_abort('invalid_input',paste('Invalid report sweep:',message))
 fields<-function(z,names_expected)if(!is.list(z)||anyDuplicated(names(z))||!setequal(names(z),names_expected))bad('missing or unexpected fields')
 if(length(x$sweeps)>20L)bad('too many sweeps')
 for(sw in x$sweeps){
  fields(sw,c('analysis','schema_version','base_spec_hash','free','rows'));i<-sw$analysis
  if(!is.numeric(i)||length(i)!=1L||!is.finite(i)||i!=floor(i)||i<1||i>length(x$specifications)||x$specifications[[i]]$analysis_kind!='power')bad('unknown power analysis')
  base<-.ig_spec_from_list(x$specifications[[i]]$specification)
  if(!identical(x$results[[i]]$status$code,'ok')||!identical(sw$schema_version,'inferguidR.report_power_sweep/0.1.0')||
   !identical(sw$base_spec_hash,base$spec_hash)||!.ig_is_string(sw$free)||!sw$free%in%names(base$inputs)||sw$free%in%base$target)bad('base or version mismatch')
  if(!is.list(sw$rows)||!is.null(names(sw$rows))||length(sw$rows)<1L||length(sw$rows)>1000L)bad('invalid row array')
  for(row in sw$rows){
   fields(row,c('value','status','reason','specification','result'))
   if(!is.numeric(row$value)||length(row$value)!=1L||!is.finite(row$value)||!.ig_is_string(row$status)||!is.character(row$reason)||length(row$reason)!=1L||is.na(row$reason))bad('invalid row metadata')
   expected<-tryCatch(stats::update(base,inputs=stats::setNames(list(row$value),sw$free)),inferguidR_invalid_input=function(err)err)
   if(inherits(expected,'condition')){
    if(!is.null(row$specification)||!is.null(row$result)||!identical(row$status,'invalid_input')||!identical(row$reason,conditionMessage(expected)))bad('invalid-input row contradicts settings')
   }else{
    spec<-.ig_spec_from_list(row$specification);result<-.ig_result_from_list(row$result);envelope<-analysis_envelope(spec,result)
    if(!identical(expected$spec_hash,envelope$spec_hash)||!identical(row$status,result$status$code)||!identical(row$reason,result$status$reason))bad('row result does not match requested settings')
   }
  }
 }
 invisible(TRUE)
}

.ig_report_sweep_figures <- function(sweeps,locale) {
 out<-list()
 for(k in seq_along(sweeps)){
  sw<-sweeps[[k]];valid<-Filter(function(row)!is.null(row$specification),sw$rows)
  if(!length(valid))next
  spec<-.ig_spec_from_list(valid[[1]]$specification);def<-.ig_get_method(spec$method_id)
  for(target in spec$target){
   data<-lapply(sw$rows,function(row){
    value<-if(identical(row$status,'ok'))row$result$solved[[target]] else NULL
    if(!is.numeric(value)||length(value)!=1L||!is.finite(value))value<-NULL
    list(x=row$value,y=value,status=row$status,reason=row$reason)
   })
   label<-function(name){x<-def$params[[name]]$label[[locale]];if(is.null(x))name else x}
   out[[length(out)+1L]]<-list(analysis=sw$analysis,plot_type='power_sweep',status=if(any(vapply(data,function(z)!is.null(z$y),TRUE)))'ok' else 'unavailable',
    sweep_index=k,target=target,axes=list(x=label(sw$free),y=label(target)),data=data,
    caption=if(locale=='de')'Einzeln berechnete Entwuerfe ohne Interpolation. Fehlgeschlagene Zeilen und Gruende bleiben im Anhang erhalten. Dies ist kein Unsicherheitsband.' else 'Individually evaluated designs without interpolation. Failed rows and reasons remain in the appendix. This is not an uncertainty band.')
  }
 }
 out
}

.ig_report_draw_sweep <- function(v,locale) {
 good<-Filter(function(z)!is.null(z$y),v$data)
 x<-vapply(good,`[[`,0,'x');y<-vapply(good,`[[`,0,'y')
 graphics::plot(x,y,pch=19,col='#304f9c',xlab=v$axes$x,ylab=v$axes$y,
  main=if(locale=='de')'Berechnete Entwurfsvarianten' else 'Evaluated design variants')
 invisible(NULL)
}
