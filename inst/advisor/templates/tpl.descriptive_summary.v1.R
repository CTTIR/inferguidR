# Synthetic summaries selected by measurement meaning, not by R storage type.
# No hypothesis test or p-value is reported.
# Continuous mean intervals assume independent observations with normal errors.
# Category intervals are marginal Clopper-Pearson intervals, not simultaneous.
analyse <- function(x, measurement, conf.level = 0.95, event = NULL) {
    spec <- infercalcR::method_spec("descriptive_summary", alternative = "two.sided", conf.level = conf.level)
    infercalcR::run_inference(spec, x = x, measurement = measurement, event = event)$payload$native
}
# With no recorded outcome type, the standalone example explicitly uses continuous data.
if(!exists("example_measurement",inherits=FALSE)) example_measurement<-"continuous"
example_x<-switch(example_measurement,
 continuous=c(1,2,4,6,7,9),
 binary=factor(c("yes","no","yes","no","yes"),levels=c("no","yes")),
 ordinal=ordered(c("low","low","mid","high"),levels=c("low","mid","high")),
 nominal=c("A","A","B","C"))
example_result<-analyse(example_x,measurement=example_measurement,conf.level=confidence_level,
 event=if(example_measurement=="binary") "yes" else NULL)
print(example_result)
