#' Export the selected advisor method as an executable R teaching example
#'
#' Returns a bundled example for the exact selected method. The script uses
#' synthetic data and records the recommendation, assumptions and rule version
#' as comments. It does not analyse user data or evaluate code during export.
#' A missing template raises an unsupported condition; no other test is used.
#'
#' The bundled examples cover one-sample, Welch, pooled two-sample and paired
#' t-tests; exact binomial and uncorrected score proportion tests; Pearson,
#' Spearman and Kendall correlation; one-sample and paired sign tests; exact
#' conditional and mid-p McNemar tests; fixed-margin Fisher tests; and
#' central, minlike and directional Boschloo tests using the package engine;
#' one-sample/paired signed-rank and independent rank-sum tests; Pearson
#' goodness-of-fit and independence tests, independent-proportion score tests,
#' and corrected/uncorrected asymptotic McNemar tests. Omnibus examples do not
#' implement directional hypotheses. Sparse-cell diagnostics never switch tests.
#' Independent-proportion intervals use unpooled normal variance, while their
#' score test uses pooled variance; decisions need not coincide. Other new
#' chi-squared examples explicitly do not compute confidence intervals.
#' General-table Fisher examples condition on both margins and order tables by
#' probability. They reject directional hypotheses, supply no sampling interval,
#' and never substitute simulation or an approximation after resource failure.
#' Their stated teaching limits are 100 units and five categories per dimension.
#' Descriptive examples require explicit continuous, binary, ordinal or nominal
#' measurement meaning, report mean or marginal category intervals and no p-values.
#' Binary inputs must declare both factor levels, or use logical values and an
#' explicit event. Ordinal quantiles are observed ordered categories (type 1).
#' Descriptive export requires alternative = NULL because it has no hypothesis.
#' Independent-group omnibus examples cover pooled and Welch one-way ANOVA and
#' Kruskal-Wallis, with explicit variance assumptions or tied-rank correction.
#' They do not supply sampling intervals or unplanned pairwise comparisons.
#' Repeated-measures examples cover tied-rank Friedman and binary Cochran Q
#' omnibus tests for three or more occasions, with confirmed unique participant IDs.
#' Missing occasions require an explicit complete-case choice and exclusions are
#' counted; validity requires justified missing-completely-at-random assumptions.
#' Zero conditional variance is unavailable. No intervals or pairwise tests are invented.
#' Paired nominal examples distinguish Stuart-Maxwell marginal homogeneity from
#' Bowker symmetry, require identical named categories on both axes and retain
#' the supplied count of excluded pairs. Singular marginal covariance or zero
#' symmetry pair totals are explicitly unavailable, without category pooling.
#' Linear regression examples use prespecified numeric predictors, an intercept
#' and classical homoscedastic t/F inference. Complete-case exclusion is explicit
#' and counted. Coefficients describe associations. Multiple-regression examples
#' separate the omnibus test from an optional prespecified nested predictor-set
#' test on exactly the same rows. No stepwise selection or robust-error switch is
#' performed. Rank-deficient and numerically exact fits cannot supply inference.
#' Model routes export two-sided inference. The standalone simple-regression
#' example also supports prespecified directional slope tests, separately from
#' the two-sided coefficient table. Correlation-to-regression overrides remain refused.
#' PCA examples form exploratory weighted sums of centered numeric variables,
#' standardised by default. The analyst chooses the retained component count;
#' weights, scores, variance fractions and reconstruction error are reported.
#' Missing-case exclusion is explicit. Repeated eigenvalues identify subspaces
#' rather than unique individual axes. Components are not latent factors.
#' PCA export requires alternative = NULL; conf.level is not used by its script.
#' Fixed-cutoff diagnostic counts examples require declared sampling, rule origin
#' and evaluation split. They retain typed undefined/infinite metrics and marginal
#' Clopper-Pearson intervals, with observed and prevalence-scenario values separate.
#' A tuned cutoff cannot claim validation on its tuning data. Four aggregate cells
#' do not identify ROC/PR curves, AUC or a cutoff sweep. Independent units and a
#' fixed rule must be justified; no clustering or selection uncertainty is inferred.
#' These examples require alternative = NULL and use conf.level for marginal intervals.
#' Brunner-Munzel examples estimate P(X > Y) plus half the tie probability,
#' using an explicit heteroscedastic placement-variance t approximation. They
#' require at least ten observations per group in this teaching implementation.
#' Zero placement variance is unavailable. Approximate intervals are intersected
#' with `[0,1]` and unconstrained limits retained. This is not a mean/median-shift test.
#' Monte Carlo goodness-of-fit examples use prespecified multinomial probabilities,
#' a recorded seed and replicate count, and restore the calling RNG state. Their
#' plus-one p-value and computational precision interval are distinct from a
#' sampling interval for an effect. Probabilities are normalised only for roundoff
#' within 64 machine epsilons; no fitted probabilities or automatic asymptotic
#' fallback are used. Simulation allocation is limited to ten million cells.
#' Cochran-Armitage trend examples require a recorded binary outcome, named
#' independent-binomial group counts and explicitly confirmed, increasing scores
#' chosen before inspecting outcomes. Score spacing is part of the hypothesis.
#' The public omnibus route exports two-sided inference; the standalone function
#' also supports a prespecified signed direction. No event variation is unavailable,
#' and sparse counts do not trigger an exact or continuity-corrected substitute.
#' No risk-difference or logistic-slope confidence interval is manufactured.
#' Jonckheere-Terpstra examples require a prespecified group order and independent
#' unit IDs, with numeric or explicitly ordered outcomes. Inclusive permutation
#' tails use a recorded seed and plus-one correction; the two-sided p-value is
#' the capped doubled smaller tail. The public route exports two-sided inference;
#' the standalone function accepts a prespecified increasing or decreasing trend.
#' Caller RNG state is restored. Simulation intervals describe computational
#' precision, not population effects. All-equal outcomes have p = 1 without
#' simulation. Limits are 200 rows, 10 groups and ten million pair comparisons.
#' Binormal classifier examples are hypothetical fixed-cutoff scenarios with
#' explicit assumed AUC, standard-deviation ratio, prevalence and count scale.
#' Higher scores are positive, on the scale standardised to non-events. Analytic
#' class-conditional rates are separate from expected-count metrics, including
#' zero-prevalence denominators. Expected N does not validate sample-size planning
#' or supply empirical confidence intervals. Finite cutoffs are required by the
#' manual-rule API. These scripts do not use conf.level and require alternative = NULL.
#' Scored classifier examples require explicit independent-unit identifiers, score direction
#' and comparator. Exploration enumerates decision sets descriptively without intervals;
#' a predeclared cutoff requires an explicit declaration. The held-out example selects
#' Youden J on tuning data with a fewest-alerts tie break, then transfers the fixed
#' numeric rule. Known overlapping IDs are refused, but disjoint IDs do not prove
#' independence. Symbolic alert-none rules without a numeric cutoff cannot transfer.
#' No score-model training or nested resampling is implied. Missing data are rejected.
#' All three require alternative = NULL; only predeclared and held-out scripts use
#' conf.level, for marginal fixed-rule Clopper-Pearson intervals.
#' ANCOVA examples fit a prespecified additive common-slope model to independent
#' units, with one categorical group and numeric covariates. The overall group F
#' test compares nested models fitted to exactly the same retained rows. Adjusted
#' means use pooled covariate means; all pairwise differences are explicitly
#' unadjusted marginal contrasts. Common slopes must be confirmed before fitting;
#' no interaction pretest or causal interpretation is implied. Rank deficiency,
#' constant responses and numerically exact fits are refused. Missing data require
#' an explicit complete-case choice. Public examples match two versus multiple
#' groups; the omnibus test remains nondirectional. Mean confidence intervals are
#' not prediction intervals for individuals.
#' Spread examples cover the normal-sample two-group F variance-ratio test,
#' median-centred Brown-Forsythe (Levene) and Fligner-Killeen tests. The recorded
#' spread route is nondirectional. The F interval concerns the raw population
#' variance ratio; robust omnibus examples do not supply sampling intervals.
#' These examples answer a spread question, never select a mean test by pretesting.
#' Their contracts reject missing or infinite
#' measurements. Rank correlations explicitly use asymptotic tests. Paired
#' examples require one wide row per confirmed, unique participant identifier.
#' Mean-test intervals concern the raw mean difference (x minus y or a stated
#' reference). Proportion intervals concern the event probability. Pearson
#' intervals use the Fisher transform; rank-correlation intervals are explicitly
#' unavailable in these base R examples. Sign and McNemar intervals concern
#' positive-sign or discordance probabilities, not median or marginal differences.
#' Mid-p examples explicitly use an exact conditional interval, not an inversion
#' of the mid-p test. Fisher intervals concern the conditional odds ratio and
#' inherit base R's numerical root tolerance. Boschloo examples preserve the
#' selected ordering and direction; sampling intervals are explicitly unavailable,
#' and numerical optimisation bounds are labelled separately. Signed-rank examples
#' require interval scale and symmetry. Their intervals concern the pseudomedian.
#' Rank-sum examples report descriptive superiority; a shift interval is available
#' only under an explicitly requested location-shift model, never for ordinal codes.
#' Exact rank examples reject zeros/ties where relevant and samples above their
#' stated teaching limit; approximation requires an explicit choice.
#' Rank calculations retain the supplied numeric precision; signed-rank
#' approximations omit null differences explicitly across R versions.
#' The pooled-Z unconditional example uses its own bounded small-sample implementation,
#' with a zero difference null, fixed group sizes and a confirmed independent binomial design.
#' Total sample size is at most 30. Two-sided inference orders absolute pooled Z.
#' Endpoints are included in the nuisance search; lower and upper numerical bounds
#' are separate from sampling uncertainty. The example refuses non-convergence and
#' withholds a decision when its numerical interval contains the significance threshold.
#' Factorial examples use two or three fixed crossed factors, full interactions,
#' sum contrasts and prespecified type-III coefficient-block F tests. Each declared
#' cell needs two observations; empty cells, repeated units and unjustified equal
#' variance are refused. Raw term p-values and Holm adjustment across all terms
#' are separate. Complete-case exclusion is explicit. Partial eta squared is not
#' an additive fraction of total explained variation in unbalanced designs.
#' Repeated one-way examples apply the Greenhouse-Geisser correction routinely
#' to both F degrees of freedom, without a Mauchly pretest. Three to thirty
#' occasions and at least three complete independent subjects are required.
#' Subject exclusion requires an explicit complete-case choice. Degenerate
#' contrast variance is refused; no post-hoc tests or effect-size interval are
#' computed. The correction is approximate, not an exact covariance-free test.
#' Regression examples use prespecified numeric predictors and explicit independent
#' units. Logistic slope odds ratios receive profile-likelihood intervals; separate
#' Wald coefficient summaries are labelled. Poisson and NB2 count models require
#' positive exposure offsets; NB2 dispersion is re-estimated in both nested models.
#' Ordered outcomes require a justified proportional-odds model, and multinomial
#' outcomes require an explicit reference category. MASS and nnet are required for
#' the corresponding examples. Singular, near-boundary and unresolved fits are
#' refused. Complete-case exclusion is explicit; no model is selected by a pretest.
#' Clustered Gaussian examples use nlme REML random-intercept models with a
#' prespecified compound-symmetric covariance and containment degrees of freedom.
#' The repeated example supports an occasion factor; the cluster example uses a
#' cluster-level factor. Binary examples distinguish conditional lme4 mixed-model
#' odds ratios from marginal geepack GEE odds ratios. GEE intervals use the robust
#' sandwich without a small-cluster correction. MAR available-case likelihood and
#' MCAR unweighted GEE assumptions are explicit. Boundary and unresolved fits are
#' refused, and the implementation cluster minimum is not an adequacy guarantee.
#' Prediction examples refit preprocessing and the prespecified model inside every
#' individual bootstrap resample. They report apparent and optimism-corrected
#' performance, Monte Carlo standard errors and all failed-resample counts.
#' Any failed resample prevents selective averaging. Logistic risk examples
#' include fixed-orientation AUC, Brier/log loss and calibration intercept/slope.
#' Internal validation is not external validation; Monte Carlo precision is not
#' a sampling interval. No cutoff, hypothesis-test power or prediction-model
#' sample-size calculation is substituted. RNG state is restored after use.
#' Common-factor examples fit a prespecified factor count by maximum likelihood
#' with oblique promax rotation and explicit Heywood-boundary refusal. They do
#' not substitute PCA or select factor counts automatically. Clustering examples
#' state k, scaling, distance and algorithm, and assess80-percent subsampling
#' stability with adjusted Rand indices. Gower/PAM requires declared variable
#' types and weights. Cluster is required for Gower/PAM. Silhouette and internal
#' stability do not establish real populations; failed subsamples are not
#' selectively averaged. Exploratory procedures have no directional alternative.
#' New-case linear predictions offer pointwise conditional-mean or new-outcome
#' Student-t intervals under independent homoscedastic normal errors. Logistic
#' risk intervals use ordinary-bootstrap percentiles (quantile type 6), requiring
#' at least five expected draws per tail; this resolution floor does not guarantee
#' coverage. Risk estimation intervals do not quantify binary-outcome uncertainty
#' or validation-performance uncertainty. All resamples must fit successfully.
#' The common-factor template also offers independent-column permutation parallel
#' analysis with SMC communalities and explicit quantile thresholds. This is
#' exploratory factor-count guidance with Monte Carlo variation, not automatic
#' model selection. Prediction-performance exports require alternative = NULL;
#' their confidence setting is used for the separate new-case prediction intervals.
#' Diagnostic plots do not trigger a switch of test. Planning availability remains the
#' availability recorded by the advisor; the script does not invent a sample
#' size or substitute a pooled test for Welch planning.
#'
#' @inheritParams advisor_recommendation
#' @param alternative `NULL` for a recorded two-sided design, or an explicit
#'   `"greater"` or `"less"` for a recorded pre-specified one-sided design.
#'   The meaning of `"greater"` is stated in the method-specific script (for
#'   example x minus y, positive-sign probability, or an odds ratio above one).
#'   It must be chosen before inspecting results.
#' @param conf.level Confidence level strictly between zero and one. This is
#'   an explicit example setting, not a value inferred from the advisor answers.
#' @return A single character string containing an R script. Use `writeLines()`
#'   to save it. Running the script prints a result; some examples also draw plots.
#' @family advisor
#' @export
advisor_r_code <- function(session, alternative = NULL, conf.level = 0.95, rules = NULL) {
  r <- advisor_recommendation(session, rules = rules)
  if (!r$status %in% c("recommended", "conditional") || is.null(r$preferred))
    .ig_abort("unsupported", "A method recommendation is required before exporting an R example.")
  path <- .ig_advisor_template_file(r$r_implementation$template_id)
  if (!nzchar(path)) .ig_abort("unsupported", "No R template is bundled for the selected variant; no substitute is generated.")
  direction <- r$design$answers[["q.direction"]]
  linear <- r$preferred$analysis_variant_id %in% c("regression_linear_simple", "regression_linear_multiple")
  # Modelling routes have no directional question; their example defaults to two-sided inference.
  if (linear && is.null(direction)) direction <- "two_sided"
  # The multi-group route records no signed direction; its trend example is two-sided.
  if (identical(r$preferred$analysis_variant_id, "jonckheere_terpstra") && is.null(direction))
    direction <- "two_sided"

  omnibus <- r$preferred$analysis_variant_id %in% c("ancova_baseline_adjusted", "chisq_goodness_of_fit", "chisq_goodness_of_fit_mc", "cochran_armitage_trend",
    "chisq_independence_rxc", "fisher_exact_rxc", "mcnemar_asymptotic_uncorrected", "mcnemar_asymptotic_cc",
    "anova_one_way", "anova_one_way_welch", "anova_factorial", "anova_repeated_one_way_gg", "regression_logistic_binary", "regression_poisson", "regression_negative_binomial", "regression_ordinal_proportional_odds", "regression_multinomial", "lmm_repeated_measures", "lmm_cluster_random_intercept", "glmm_logistic_random_intercept", "gee_logistic_exchangeable", "kruskal_wallis", "friedman", "cochran_q",
    "stuart_maxwell_marginal_homogeneity", "bowker_symmetry")
  prediction <- r$preferred$analysis_variant_id %in% c("prediction_model_linear", "prediction_model_logistic")
  spread <- r$preferred$analysis_variant_id %in% c("variance_f_test_two",
    "variance_brown_forsythe", "variance_fligner_killeen")
  descriptive <- identical(r$preferred$analysis_variant_id, "descriptive_summary")
  exploratory <- r$preferred$analysis_variant_id %in% c("pca_components", "efa_common_factor_ml", "cluster_hierarchical_ward", "cluster_kmeans", "cluster_gower_pam")
  fixed_counts <- identical(r$preferred$analysis_variant_id, "diagnostic_2x2_fixed_cutoff")
  theoretical <- identical(r$preferred$analysis_variant_id, "classifier_binormal_theory")
  scored_classifier <- r$preferred$analysis_variant_id %in% c("classifier_cutoff_exploration", "classifier_predeclared_rule", "classifier_heldout_validation")
  exploration <- identical(r$preferred$analysis_variant_id, "classifier_cutoff_exploration")
  if (prediction) {
    if (!is.null(alternative)) .ig_abort("invalid_input", "Prediction performance has no hypothesis alternative.")
  } else if (scored_classifier) {
    if (!is.null(alternative)) .ig_abort("invalid_input", "Classifier operating metrics have no hypothesis alternative.")
  } else if (theoretical) {
    if (!is.null(alternative)) .ig_abort("invalid_input", "A theoretical cutoff scenario has no hypothesis alternative.")
  } else if (fixed_counts) {
    if (!is.null(alternative)) .ig_abort("invalid_input", "Diagnostic operating metrics have no hypothesis alternative.")
  } else if (exploratory) {
    if (!is.null(alternative)) .ig_abort("invalid_input", "The selected exploratory procedure has no hypothesis alternative.")
  } else if (descriptive) {
    if (!is.null(alternative))
      .ig_abort("invalid_input", "Descriptive summaries have no hypothesis alternative.")
  } else if (spread) {
    if (!is.null(direction) && !identical(direction, "two_sided") ||
        (!is.null(alternative) && !identical(alternative, "two.sided")))
      .ig_abort("invalid_input", "The recorded spread-comparison route supports only a nondirectional comparison.")
    alternative <- "two.sided"
  } else if (omnibus) {
    if (identical(direction, "one_sided_prespecified") ||
        (!is.null(alternative) && !identical(alternative, "two.sided")))
      .ig_abort("invalid_input", "The selected example is omnibus; it cannot implement a directional hypothesis.")
    alternative <- "two.sided"
  } else if (identical(direction, "two_sided")) {
    if (is.null(alternative)) alternative <- "two.sided"
    if (!identical(alternative, "two.sided"))
      .ig_abort("invalid_input", "The requested alternative conflicts with the recorded two-sided design.")
  } else if (identical(direction, "one_sided_prespecified")) {
    if (!.ig_is_string(alternative) || !alternative %in% c("greater", "less"))
      .ig_abort("invalid_input", "Specify greater or less from the pre-specified hypothesis; the advisor does not record its sign.")
  } else .ig_abort("invalid_input", "Resolve the hypothesis direction in the advisor before exporting code.")
  measurement_setting <- NULL
  if (identical(r$preferred$analysis_variant_id, "ancova_baseline_adjusted"))
    measurement_setting <- paste0("example_group_count <- ", if(identical(r$design$answers[["q.n_groups"]], "more_than_two")) "3L" else "2L")
  if (descriptive) {
    measurement <- r$design$answers[["q.outcome_type"]] %||% "continuous"
    if (!measurement %in% c("continuous", "binary", "ordinal", "nominal"))
      .ig_abort("unsupported", "No descriptive example is bundled for the recorded measurement type.")
    measurement_setting <- paste0("example_measurement <- ", encodeString(measurement, quote='"'))
  }
  if (!is.numeric(conf.level) || length(conf.level) != 1L || !is.finite(conf.level) || conf.level <= 0 || conf.level >= 1)
    .ig_abort("invalid_input", "conf.level must be a finite number strictly between zero and one.")
  comment <- function(value) paste0("# ", gsub("[[:cntrl:]]", " ", value))
  lang <- session$language
  assumption_text <- vapply(r$assumptions, function(a) paste0(
    "[", a$state, "; ", a$source, "; ", a$stated_for %||% "design", "] ",
    a[[paste0("text_",lang)]] %||% a$text_en), "")
  body <- readLines(path, warn = FALSE, encoding = "UTF-8")
  header <- c("inferguidR synthetic R analysis example", r$preferred[[paste0("name_",lang)]],
    paste("Variant:", r$preferred$analysis_variant_id), paste("Recommendation:", r$recommendation_hash),
    paste("Status:", r$status), paste("Rules:", r$rules$rules_version),
    paste("Template SHA-256:", digest::digest(file=path,algo="sha256")),
    "Recorded design and assumptions have not been checked against data.",
    r$rationale[[lang]] %||% r$rationale$en,
    paste("Estimand:", r$estimand[[lang]] %||% r$estimand$en),
    vapply(names(r$design$answers), function(id) paste(id, r$design$answers[[id]], sep = " = "), ""),
    assumption_text, paste("Planning:", r$planning[[paste0("note_",lang)]] %||% r$planning$note_en))
  paste(c(vapply(header, comment, ""),
    measurement_setting,
    if (!descriptive && !exploratory && !fixed_counts && !theoretical && !scored_classifier && !prediction) paste0("alternative <- ", encodeString(alternative, quote='"')),
    if (!exploratory && !theoretical && !exploration) paste0("confidence_level <- ", format(conf.level, digits = 17, scientific = TRUE, decimal.mark = ".")), body), collapse="\n")
}
