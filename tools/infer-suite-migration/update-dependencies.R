# Refresh the direct dependency inventory from the three actual declarations.
# Execution outcomes belong to a dated validation record, not an installation probe.
args <- commandArgs(TRUE)
stopifnot(length(args) == 3L)
packages <- c(infercalcR=args[1], inferguidR=args[2], infereportR=args[3])
base <- c('base','compiler','datasets','grDevices','graphics','grid','methods','parallel','splines','stats','stats4','tcltk','tools','utils')
recommended <- c('boot','class','cluster','codetools','foreign','KernSmooth','lattice','MASS','Matrix','mgcv','nlme','nnet','rpart','spatial','survival')
references <- c('Exact','exact2x2','MESS','pROC','pwr','pwrss','WebPower','yardstick')
feature <- c(infercalcR='Shared numerical provider and result contract',
  digest='Historical SHA256 identity and source data fingerprint',
  jsonlite='Historical canonical JSON and correctly rounded decimal parsing',
  cluster='Gower distance and PAM clustering', nlme='Gaussian mixed models',
  MASS='Negative-binomial and ordinal regression', nnet='Multinomial regression',
  lme4='Logistic random-intercept mixed models', geepack='Exchangeable logistic GEE',
  ggplot2='Finalized-result graphical annotations',
  DescrTab2='Pinned precomputed native table import and render',
  jsonvalidate='Optional historical JSON schema validation',
  shiny='Guidance workbench application', bslib='Guidance workbench presentation',
  knitr='Vignette construction', rmarkdown='Vignette construction',
  testthat='Retained legacy regression coverage', withr='Retained test environment isolation',
  shinytest2='Existing workbench browser regression', quarto='Existing optional document tooling')
rows <- list()
for (pkg in names(packages)) {
  d <- read.dcf(file.path(packages[[pkg]],'DESCRIPTION'))
  for (field in intersect(c('Imports','Suggests'),colnames(d))) {
    deps <- trimws(gsub(' *\\([^)]*\\)','',strsplit(d[1,field],',')[[1]]))
    for (dep in deps) {
      category <- if (dep %in% base) 'base' else if (dep %in% recommended) 'recommended' else 'contributed'
      role <- if (dep %in% references) 'optional independent reference comparator' else if (dep %in% c('testthat','withr','shinytest2')) 'development test' else if (dep %in% c('knitr','rmarkdown')) 'development documentation' else 'runtime'
      why <- if (dep %in% names(feature)) unname(feature[dep]) else if (dep %in% references) 'Separate historical reference comparison, no production method substitution' else 'R-distributed numerical, graphical or utility functionality'
      version <- tryCatch(as.character(utils::packageVersion(dep)), error=function(e) 'NOT_INSTALLED')
      rows[[length(rows)+1L]] <- data.frame(package=pkg,dependency=dep,category=category,role=role,mandatory=field=='Imports',
        feature=why,base_R_alternative=if(category=='base') 'standard R' else 'No silent substitution of method or historical identity',
        disposition=if(field=='Imports') 'Imports' else 'Suggests with feature-local use',justification=why,
        version_tested=version,validation_state='See dated installed-candidate evidence',stringsAsFactors=FALSE)
    }
  }
}
write.csv(do.call(rbind,rows), 'tools/infer-suite-migration/dependencies.csv', row.names=FALSE)
