# Synthetic independent-binomial example: row1 event probability minus row2.
# Column1 is event and column2 non-event. Fix direction before observing data.
# Uncorrected pooled-score test; the interval uses unpooled estimated variance.
analyse <- function(counts, alternative = "two.sided", conf.level = 0.95) {
    spec <- infercalcR::method_spec("prop_two_independent_pearson_chisq", alternative = alternative,
        conf.level = conf.level)
    infercalcR::run_inference(spec, counts = counts)$payload$native
}
example_counts<-rbind(Group1=c(Event=30,Non_event=70),Group2=c(Event=45,Non_event=55))
example_result<-analyse(example_counts,alternative,confidence_level)
print(example_result)
