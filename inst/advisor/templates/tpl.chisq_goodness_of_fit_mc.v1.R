# Independent multinomial counts against probabilities specified before seeing the data.
# Monte Carlo uncertainty is computational, not a sampling interval for population effects.
analyse <- function(observed, probabilities, seed, B = 19999L, conf.level = 0.95, mc.level = 0.99) {
    spec <- infercalcR::method_spec("chisq_goodness_of_fit_mc", alternative = "two.sided", conf.level = conf.level)
    infercalcR::run_inference(spec, observed = observed, probabilities = probabilities, seed = seed,
        B = B, mc.level = mc.level)$payload$native
}
example_result<-analyse(c(A=1,B=3,C=6),c(A=.2,B=.3,C=.5),seed=104729L,B=19999L,conf.level=confidence_level)
print(example_result)
