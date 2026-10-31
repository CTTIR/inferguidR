args<-commandArgs(TRUE);stopifnot(length(args)==3L)
p<-normalizePath(args[1],mustWork=TRUE);out<-normalizePath(args[2],mustWork=TRUE);lane<-args[3]
options(Ncpus=1)
library(inferguidR)
cat('Installed guid:',find.package('inferguidR'),'\nInstalled calc:',find.package('infercalcR'),'\n')
writeLines(capture.output(sessionInfo()),file.path(out,paste0(lane,'-session.txt')))
filter<-switch(lane,
 core='^(core-|json-number-precision|boschloo-api$|boschloo-design$|boschloo-method$|classif-spec$|classif-release$|classif-export$|num-dist$)',
 legacy='^(app-|boschloo-enumeration$|boschloo-exact$|boschloo-plotdata$|boschloo-properties$|boschloo-report$|classif-|compromise-|glm-|gof-|learning-|methods-|num-large|num-modes|num-root|power-|protocol-|registry-|report-|slice-|sweep-|unconditional-|variance-|workbench-)',
 stop('unknown lane'))
# Work from outside a source tree and explicitly load the installed namespace.
setwd(tempdir())
r<-testthat::test_dir(p,package='inferguidR',load_package='installed',filter=filter,reporter='summary',stop_on_failure=FALSE)
d<-as.data.frame(r)
write.csv(d[,setdiff(names(d),'result')],file.path(out,paste0(lane,'-tests.csv')),row.names=FALSE)
stopifnot(nrow(d)>0L,!any(d$error),sum(d$failed)==0L)
cat('PASS',nrow(d),'blocks',sum(d$passed),'assertions;',sum(d$skipped),'skips;',sum(d$warning),'warnings.\n')
