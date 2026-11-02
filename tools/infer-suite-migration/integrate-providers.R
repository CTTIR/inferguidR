# Development-time migration from the recorded symbol map; never used at runtime.
args<-commandArgs(TRUE);stopifnot(length(args)==1L)
calc<-normalizePath(args[1],mustWork=TRUE)
map<-read.csv(file.path(calc,'tools/extraction/source-map.csv'),stringsAsFactors=FALSE)
# Historical schema validation and local condition construction remain guidance-owned.
structural<-map$symbol[map$source_file=='core-envelope.R']
local_status<-map$symbol[map$source_file=='core-status.R' & !map$symbol%in%c('ig_status','ig_status_codes','is_ig_status','format.ig_status','print.ig_status')]
local_simple<-c('.bl_is_count','.bl_check_alpha','.bl_check_prob','.bl_fmt')
keep<-c(structural,local_status,local_simple)
remove<-setdiff(map$symbol,keep)
removed<-list()
files <- unique(paste0('R/',map$source_file))
for(f in files) {
 if(basename(f)=='suite-results.R')next
 txt<-system2('git',c('show',paste0('aad7fd6:',f)),stdout=TRUE)
 if(!length(txt))next
 p<-parse(text=txt,keep.source=TRUE);drop<-integer();nms<-character()
 for(i in seq_along(p)){
  e<-p[[i]]
  if(!is.call(e)||!as.character(e[[1]])%in%c('<-','='))next
  lhs<-e[[2]];nm<-if(is.symbol(lhs)||is.character(lhs))as.character(lhs)else if(is.call(lhs)&&identical(lhs[[1]],as.name('$')))as.character(lhs[[2]])else ''
  if(nm%in%remove){
   sr<-attr(p,'srcref')[[i]]
   first<-if(i==1L)1L else attr(p,'srcref')[[i-1L]][3]+1L
   last<-if(i==length(p))length(txt)else sr[3]
   drop<-union(drop,seq.int(first,last));nms<-c(nms,nm)
  }
 }
 if(!length(drop))next
 left<-txt[-drop]
 # Compact the whitespace left by the deleted definitions and their comments.
 left<-strsplit(gsub('\n[ \t]*\n[ \t]*\n+','\n\n',paste(left,collapse='\n')),'\n',fixed=TRUE)[[1]]
 if(!length(parse(text=left)))unlink(f) else writeLines(left,f)
 removed[[f]]<-data.frame(source_file=f,symbol=nms,treatment='public calc provider',stringsAsFactors=FALSE)
}
forward<-readLines(file.path(calc,'tools/extraction/guid-forwarders.R'))
# Keep the historical schema integrity validator with its guid reader; hashes delegate.
start<-grep('^\\.ig_check_aux_spec <-',forward)
if(length(start))forward<-forward[-start]
# Description of actual runtime architecture, not implementation instructions.
forward[1:2]<-c('# Guidance consumes canonical numerical and identity providers from infercalcR.', '# Historical schema dispatch remains local; no numerical engines are duplicated.')
writeLines(forward,'R/suite-providers.R')
ns<-readLines('NAMESPACE');oldexp<-sub('export\\((.*)\\)','\\1',ns[grepl('^export\\(',ns)])
movedexp<-intersect(oldexp,remove)
calcns<-readLines(file.path(calc,'NAMESPACE'));s3<-calcns[grepl('^S3method',calcns)]
ns<-ns[!ns%in%s3]
ns<-unique(c(ns,paste0('importFrom(infercalcR,',movedexp,')')))
writeLines(ns,'tools/infer-suite-migration/provider-integration/NAMESPACE.proposed')
writeLines(paste0('importFrom(infercalcR,',movedexp,')'),'tools/infer-suite-migration/provider-integration/NAMESPACE.imports')
reexports<-c("#' Calculation API reexports", "#'", "#' Preserve the historical guidance help pages while registering their", "#' canonical calculation implementations from infercalcR.", "#' @noRd", unlist(lapply(movedexp,function(n)c(paste0("#' @rawNamespace importFrom(infercalcR,",n,")"),paste0("#' @rawNamespace export(",n,")")))), 'NULL')
writeLines(reexports,'R/reexports.R')
writeLines(s3[s3%in%readLines('NAMESPACE')],'tools/infer-suite-migration/provider-integration/NAMESPACE.remove-S3')
write.csv(do.call(rbind,removed),'tools/infer-suite-migration/provider-integration/removed-symbols.csv',row.names=FALSE)
write.csv(map[map$symbol%in%keep,], 'tools/infer-suite-migration/provider-integration/retained-structural-symbols.csv',row.names=FALSE)
cat('Removed',length(remove),'mapped objects; re-exported',length(movedexp),'legacy APIs through explicit imports.\n')
