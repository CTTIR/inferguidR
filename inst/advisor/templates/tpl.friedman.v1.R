# Synthetic complete-block teaching example, with one wide row per confirmed participant.
# Columns are prespecified occasions; participant IDs establish matching explicitly.
# Missing occasions require an explicit complete-case choice; excluded blocks are counted.
analyse <- function(data, alternative = "two.sided", conf.level = 0.95, missing = "error") {
    spec <- infercalcR::method_spec("friedman", alternative = alternative, conf.level = conf.level)
    infercalcR::run_inference(spec, data = data, missing = missing)$payload$native
}
example_data<-data.frame(id=paste0("participant_",1:6),a=c(1,2,2,4,3,2), b=c(2,1,3,4,2,3), c=c(3,3,1,4,1,1))
example_result<-analyse(example_data,alternative,confidence_level)
print(example_result)
