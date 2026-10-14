# Curated G*Power fixture comparisons for group G10 (Wilcoxon tests, fast subset) ------------
#
# A fast, offline subset of the emulator-capture comparisons of
# validation/compare/<method_id>/ for t_wilcoxon_signed_rank_one_sample,
# t_wilcoxon_signed_rank_matched_pairs and t_wilcoxon_mann_whitney (the full
# runs: validation/tools/compare_fixtures.R). Every case is a fixture record of
# validation/fixtures/emulator_capture/ or emulator_capture_supplement/g10/,
# cited by its record id; the inputs are the specification the comparator's
# mapping builds from the record and the displayed strings are copies of the
# record (G*Power 3.1.9.7 under Wine, evidence kind emulator_capture,
# decisions D22/D23). The last test checks the strings against the fixture
# files in a development checkout. One case per method x option branch x mode
# was taken from the records that agree in the full comparison; the records
# that disagree by design (inconsistent manual moments, the a priori N rule,
# the sign of d in the p1 branch) or through a reference defect or a pending
# core change are documented in validation/compare/<method_id>/summary.md and
# tested separately in test-methods-g10-wilcoxon.R.
#
# Match rules (validation/tolerances.md, fixed before comparison): integers
# exact (MR-EXACT); displayed reals within the closed rounding interval of the
# last displayed digit plus 8 * eps * max(1, |x|) (MR-DISPLAY). G*Power
# refusals ('Input error') must be matched by a typed refusal without numbers
# (TD-0010). Display conventions as in the mapping (wilcoxon_mapping_common.R
# M3-M7, M12).

ref_reset_registry()

g10_fixture_cases <- list(
  "GPW-T-WSR1-APR-B001-O01--r01" = list(method = "t_wilcoxon_signed_rank_one_sample", mode = "a_priori", branch = "A.R.E. parent",
    inputs = list(d = 0.5, alpha = 0.05, power = 0.8),
    options = list(method = "are", are_specification = "predefined_parent", tails = "one", parent_distribution = "normal"),
    outputs = c("Noncentrality parameter \u03b4" = "2.5854415", "Critical t" = "1.7062592", "Df" = "25.7380304", "Total sample size" = "28", "Actual power" = "0.8083058")),
  "GPW-T-WSR1-APR-B005-D01--r01" = list(method = "t_wilcoxon_signed_rank_one_sample", mode = "a_priori", branch = "A.R.E. manual",
    inputs = list(are = 0.2, d = 0.02, alpha = 0.001, power = 0.99),
    options = list(method = "are", are_specification = "manual_are", tails = "one"),
    outputs = c("Noncentrality parameter \u03b4" = "?", "Critical t" = "?", "Df" = "?", "Total sample size" = "?", "Actual power" = "?"),
    message = "The minimal A.R.E. of the Wilcoxon test versus the t test is 0.864."),
  "GPW-T-WSR1-APR-B009-O01--r01" = list(method = "t_wilcoxon_signed_rank_one_sample", mode = "a_priori", branch = "Lehmann d",
    inputs = list(d = 0.5, alpha = 0.05, power = 0.8),
    options = list(method = "lehmann", lehmann_input_method = "predefined_parent", effect_size_parameter = "d", continuity_correction = "off", distribution_plot = "normalized_z", tails = "one", parent_distribution = "normal"),
    outputs = c("Critical z" = "1.6448536", "Total sample size" = "28", "Actual power" = "0.8078139", "Moment p1" = "0.6914625", "Moment p2" = "0.7602499", "Moment p3" = "0.6337020")),
  "SUP-G10-T-WSR1-APR-S02--r01" = list(method = "t_wilcoxon_signed_rank_one_sample", mode = "a_priori", branch = "Lehmann moments +cc",
    inputs = list(p3 = 0.68, p2 = 0.8, p1 = 0.7, alpha = 0.05, power = 0.9),
    options = list(method = "lehmann", lehmann_input_method = "manual_moments", continuity_correction = "on", distribution_plot = "normalized_z", tails = "two"),
    outputs = c("Critical z" = "1.9599640", "Total sample size" = "32", "Actual power" = "0.9052746")),
  "GPW-T-WSR1-APR-B025-O01--r01" = list(method = "t_wilcoxon_signed_rank_one_sample", mode = "a_priori", branch = "Lehmann moments",
    inputs = list(p3 = 0.8206572, p2 = 0.8890997, p1 = 0.8387046, alpha = 0.05, power = 0.8),
    options = list(method = "lehmann", lehmann_input_method = "manual_moments", continuity_correction = "off", distribution_plot = "normalized_z", tails = "one"),
    outputs = c("Critical z" = "1.6448536", "Total sample size" = "11", "Actual power" = "0.8308343")),
  "GPW-T-WSR1-APR-B053-O01--r01" = list(method = "t_wilcoxon_signed_rank_one_sample", mode = "a_priori", branch = "Lehmann p1 +cc",
    inputs = list(relative_sd = 1, p1 = 0.6382, alpha = 0.05, power = 0.8),
    options = list(method = "lehmann", lehmann_input_method = "predefined_parent", effect_size_parameter = "p1", continuity_correction = "on", distribution_plot = "normalized_z", tails = "one", parent_distribution = "normal"),
    outputs = c("Critical z" = "1.6448536", "Total sample size" = "54", "Actual power" = "0.8037432", "Effect size d" = "0.3536516", "Moment p2" = "0.6915114", "Moment p3" = "0.5463044")),
  "GPW-T-WSR1-APR-B054-O01--r01" = list(method = "t_wilcoxon_signed_rank_one_sample", mode = "a_priori", branch = "Lehmann p1",
    inputs = list(relative_sd = 1, p1 = 0.6382, alpha = 0.05, power = 0.8),
    options = list(method = "lehmann", lehmann_input_method = "predefined_parent", effect_size_parameter = "p1", continuity_correction = "off", distribution_plot = "normalized_z", tails = "one", parent_distribution = "normal"),
    outputs = c("Critical z" = "1.6448536", "Total sample size" = "54", "Actual power" = "0.8050557", "Effect size d" = "0.3536516", "Moment p2" = "0.6915114", "Moment p3" = "0.5463044")),
  "GPW-T-WSR1-APR-B065-O01--r01" = list(method = "t_wilcoxon_signed_rank_one_sample", mode = "a_priori", branch = "Lehmann d +cc /T",
    inputs = list(d = 0.5, alpha = 0.05, power = 0.8),
    options = list(method = "lehmann", lehmann_input_method = "predefined_parent", effect_size_parameter = "d", continuity_correction = "on", distribution_plot = "normal_rank_sum", tails = "one", parent_distribution = "normal"),
    outputs = c("Lower critical T" = "275.7332", "Upper critical T" = "275.7332", "Total sample size" = "28", "Actual power" = "0.8040183", "Moment p1" = "0.6914625", "Moment p2" = "0.7602499", "Moment p3" = "0.6337020")),
  "GPW-T-WSR1-APR-B083-O01--r01" = list(method = "t_wilcoxon_signed_rank_one_sample", mode = "a_priori", branch = "Lehmann moments +cc /T",
    inputs = list(p3 = 0.8206572, p2 = 0.8890997, p1 = 0.8387046, alpha = 0.05, power = 0.8),
    options = list(method = "lehmann", lehmann_input_method = "manual_moments", continuity_correction = "on", distribution_plot = "normal_rank_sum", tails = "one"),
    outputs = c("Lower critical T" = "52.0000337", "Upper critical T" = "52.0000337", "Total sample size" = "11", "Actual power" = "0.8119807")),
  "GPW-T-WSR1-COM-B001-O01--r01" = list(method = "t_wilcoxon_signed_rank_one_sample", mode = "compromise", branch = "A.R.E. parent",
    inputs = list(d = 0.5, error_ratio = 1, n_total = 50),
    options = list(method = "are", are_specification = "predefined_parent", tails = "one", parent_distribution = "normal"),
    outputs = c("Noncentrality parameter \u03b4" = "3.4549415", "Critical t" = "1.7365483", "Df" = "46.7464829", "\u03b1 err prob" = "0.0445268", "\u03b2 err prob" = "0.0445268", "Power (1-\u03b2 err prob)" = "0.9554732")),
  "GPW-T-WSR1-COM-B005-D01--r01" = list(method = "t_wilcoxon_signed_rank_one_sample", mode = "compromise", branch = "A.R.E. manual",
    inputs = list(are = 0.2, d = 0.02, error_ratio = 0.01, n_total = 100000),
    options = list(method = "are", are_specification = "manual_are", tails = "one"),
    outputs = c("Noncentrality parameter \u03b4" = "?", "Critical t" = "?", "Df" = "?", "\u03b1 err prob" = "?", "\u03b2 err prob" = "?", "Power (1-\u03b2 err prob)" = "?"),
    message = "The minimal A.R.E. of the Wilcoxon test versus the t test is 0.864."),
  "GPW-T-WSR1-COM-B009-O01--r01" = list(method = "t_wilcoxon_signed_rank_one_sample", mode = "compromise", branch = "Lehmann d",
    inputs = list(d = 0.5, error_ratio = 1, n_total = 50),
    options = list(method = "lehmann", lehmann_input_method = "predefined_parent", effect_size_parameter = "d", continuity_correction = "off", distribution_plot = "normalized_z", tails = "one", parent_distribution = "normal"),
    outputs = c("Critical z" = "1.7401320", "\u03b1 err prob" = "0.0409179", "\u03b2 err prob" = "0.0409179", "Power (1-\u03b2 err prob)" = "0.9590821", "Moment p1" = "0.6914625", "Moment p2" = "0.7602499", "Moment p3" = "0.6337020")),
  "GPW-T-WSR1-COM-B025-O01--r01" = list(method = "t_wilcoxon_signed_rank_one_sample", mode = "compromise", branch = "Lehmann moments",
    inputs = list(p3 = 0.8206572, p2 = 0.8890997, p1 = 0.8387046, error_ratio = 1, n_total = 50),
    options = list(method = "lehmann", lehmann_input_method = "manual_moments", continuity_correction = "off", distribution_plot = "normalized_z", tails = "one"),
    outputs = c("Critical z" = "2.9694897", "\u03b1 err prob" = "0.001491474", "\u03b2 err prob" = "0.001491474", "Power (1-\u03b2 err prob)" = "0.9985085")),
  "GPW-T-WSR1-COM-B053-O01--r01" = list(method = "t_wilcoxon_signed_rank_one_sample", mode = "compromise", branch = "Lehmann p1 +cc",
    inputs = list(relative_sd = 1, p1 = 0.6382, error_ratio = 1, n_total = 50),
    options = list(method = "lehmann", lehmann_input_method = "predefined_parent", effect_size_parameter = "p1", continuity_correction = "on", distribution_plot = "normalized_z", tails = "one", parent_distribution = "normal"),
    outputs = c("Critical z" = "1.2205015", "\u03b1 err prob" = "0.1111374", "\u03b2 err prob" = "0.1111374", "Power (1-\u03b2 err prob)" = "0.8888626", "Effect size d" = "0.3536516", "Moment p2" = "0.6915114", "Moment p3" = "0.5463044")),
  "GPW-T-WSR1-COM-B054-O01--r01" = list(method = "t_wilcoxon_signed_rank_one_sample", mode = "compromise", branch = "Lehmann p1",
    inputs = list(relative_sd = 1, p1 = 0.6382, error_ratio = 1, n_total = 50),
    options = list(method = "lehmann", lehmann_input_method = "predefined_parent", effect_size_parameter = "p1", continuity_correction = "off", distribution_plot = "normalized_z", tails = "one", parent_distribution = "normal"),
    outputs = c("Critical z" = "1.2230335", "\u03b1 err prob" = "0.1106585", "\u03b2 err prob" = "0.1106585", "Power (1-\u03b2 err prob)" = "0.8893415", "Effect size d" = "0.3536516", "Moment p2" = "0.6915114", "Moment p3" = "0.5463044")),
  "GPW-T-WSR1-COM-B065-O01--r01" = list(method = "t_wilcoxon_signed_rank_one_sample", mode = "compromise", branch = "Lehmann d +cc /T",
    inputs = list(d = 0.5, error_ratio = 1, n_total = 50),
    options = list(method = "lehmann", lehmann_input_method = "predefined_parent", effect_size_parameter = "d", continuity_correction = "on", distribution_plot = "normal_rank_sum", tails = "one", parent_distribution = "normal"),
    outputs = c("Lower critical T" = "817.9888", "Upper critical T" = "817.9888", "\u03b1 err prob" = "0.0411510", "\u03b2 err prob" = "0.0411510", "Power (1-\u03b2 err prob)" = "0.9588490", "Moment p1" = "0.6914625", "Moment p2" = "0.7602499", "Moment p3" = "0.6337020")),
  "GPW-T-WSR1-COM-B083-O01--r01" = list(method = "t_wilcoxon_signed_rank_one_sample", mode = "compromise", branch = "Lehmann moments +cc /T",
    inputs = list(p3 = 0.8206572, p2 = 0.8890997, p1 = 0.8387046, error_ratio = 1, n_total = 50),
    options = list(method = "lehmann", lehmann_input_method = "manual_moments", continuity_correction = "on", distribution_plot = "normal_rank_sum", tails = "one"),
    outputs = c("Lower critical T" = "945.303", "Upper critical T" = "945.303", "\u03b1 err prob" = "0.001506142", "\u03b2 err prob" = "0.001506142", "Power (1-\u03b2 err prob)" = "0.9984939")),
  "GPW-T-WSR1-CRI-B001-O01--r01" = list(method = "t_wilcoxon_signed_rank_one_sample", mode = "criterion", branch = "A.R.E. parent",
    inputs = list(d = 0.5, power = 0.8, n_total = 50),
    options = list(method = "are", are_specification = "predefined_parent", tails = "one", parent_distribution = "normal"),
    outputs = c("Noncentrality parameter \u03b4" = "3.4549415", "Critical t" = "2.5975827", "Df" = "46.7464829", "\u03b1 err prob" = "0.006257226")),
  "GPW-T-WSR1-CRI-B005-D01--r01" = list(method = "t_wilcoxon_signed_rank_one_sample", mode = "criterion", branch = "A.R.E. manual",
    inputs = list(are = 0.2, d = 0.02, power = 0.9999, n_total = 100000),
    options = list(method = "are", are_specification = "manual_are", tails = "one"),
    outputs = c("Noncentrality parameter \u03b4" = "?", "Critical t" = "?", "Df" = "?", "\u03b1 err prob" = "?"),
    message = "The minimal A.R.E. of the Wilcoxon test versus the t test is 0.864."),
  "GPW-T-WSR1-CRI-B015-O01--r01" = list(method = "t_wilcoxon_signed_rank_one_sample", mode = "criterion", branch = "Lehmann d",
    inputs = list(d = 0.5, power = 0.8, n_total = 50),
    options = list(method = "lehmann", lehmann_input_method = "predefined_parent", effect_size_parameter = "d", continuity_correction = "off", distribution_plot = "normalized_z", tails = "one", parent_distribution = "logistic"),
    outputs = c("Critical z" = "2.6432172", "\u03b1 err prob" = "0.004106117", "Moment p1" = "0.7123653", "Moment p2" = "0.7726573", "Moment p3" = "0.6515369")),
  "GPW-T-WSR1-CRI-B025-O01--r01" = list(method = "t_wilcoxon_signed_rank_one_sample", mode = "criterion", branch = "Lehmann moments",
    inputs = list(p3 = 0.8206572, p2 = 0.8890997, p1 = 0.8387046, power = 0.8, n_total = 50),
    options = list(method = "lehmann", lehmann_input_method = "manual_moments", continuity_correction = "off", distribution_plot = "normalized_z", tails = "one"),
    outputs = c("Critical z" = "4.2558879", "\u03b1 err prob" = "1.041105e-05")),
  "GPW-T-WSR1-CRI-B053-O01--r01" = list(method = "t_wilcoxon_signed_rank_one_sample", mode = "criterion", branch = "Lehmann p1 +cc",
    inputs = list(relative_sd = 1, p1 = 0.6382, power = 0.8, n_total = 50),
    options = list(method = "lehmann", lehmann_input_method = "predefined_parent", effect_size_parameter = "p1", continuity_correction = "on", distribution_plot = "normalized_z", tails = "one", parent_distribution = "normal"),
    outputs = c("Critical z" = "1.5638522", "\u03b1 err prob" = "0.0589261", "Effect size d" = "0.3536516", "Moment p2" = "0.6915114", "Moment p3" = "0.5463044")),
  "GPW-T-WSR1-CRI-B054-O01--r01" = list(method = "t_wilcoxon_signed_rank_one_sample", mode = "criterion", branch = "Lehmann p1",
    inputs = list(relative_sd = 1, p1 = 0.6382, power = 0.8, n_total = 50),
    options = list(method = "lehmann", lehmann_input_method = "predefined_parent", effect_size_parameter = "p1", continuity_correction = "off", distribution_plot = "normalized_z", tails = "one", parent_distribution = "normal"),
    outputs = c("Critical z" = "1.5686788", "\u03b1 err prob" = "0.0583614", "Effect size d" = "0.3536516", "Moment p2" = "0.6915114", "Moment p3" = "0.5463044")),
  "GPW-T-WSR1-CRI-B065-O01--r01" = list(method = "t_wilcoxon_signed_rank_one_sample", mode = "criterion", branch = "Lehmann d +cc /T",
    inputs = list(d = 0.5, power = 0.8, n_total = 50),
    options = list(method = "lehmann", lehmann_input_method = "predefined_parent", effect_size_parameter = "d", continuity_correction = "on", distribution_plot = "normal_rank_sum", tails = "one", parent_distribution = "normal"),
    outputs = c("Lower critical T" = "894.2425", "Upper critical T" = "894.2425", "\u03b1 err prob" = "0.006688345", "Moment p1" = "0.6914625", "Moment p2" = "0.7602499", "Moment p3" = "0.6337020")),
  "GPW-T-WSR1-CRI-B083-O01--r01" = list(method = "t_wilcoxon_signed_rank_one_sample", mode = "criterion", branch = "Lehmann moments +cc /T",
    inputs = list(p3 = 0.8206572, p2 = 0.8890997, p1 = 0.8387046, power = 0.8, n_total = 50),
    options = list(method = "lehmann", lehmann_input_method = "manual_moments", continuity_correction = "on", distribution_plot = "normal_rank_sum", tails = "one"),
    outputs = c("Lower critical T" = "1078.375", "Upper critical T" = "1078.375", "\u03b1 err prob" = "1.0638e-05")),
  "SUP-G10-T-WSR1-POS-D01--r01" = list(method = "t_wilcoxon_signed_rank_one_sample", mode = "post_hoc", branch = "A.R.E. parent",
    inputs = list(d = 0.5, alpha = 0.05, n_total = 3),
    options = list(method = "are", are_specification = "predefined_parent", tails = "one", parent_distribution = "normal"),
    outputs = c("Noncentrality parameter \u03b4" = "0.8462844", "Critical t" = "3.0696211", "Df" = "1.8647890", "Power (1-\u03b2 err prob)" = "0.1450539")),
  "GPW-T-WSR1-POS-B005-D01--r01" = list(method = "t_wilcoxon_signed_rank_one_sample", mode = "post_hoc", branch = "A.R.E. manual",
    inputs = list(are = 0.2, d = 0.02, alpha = 1e-04, n_total = 100000),
    options = list(method = "are", are_specification = "manual_are", tails = "one"),
    outputs = c("Noncentrality parameter \u03b4" = "?", "Critical t" = "?", "Df" = "?", "Power (1-\u03b2 err prob)" = "?"),
    message = "The minimal A.R.E. of the Wilcoxon test versus the t test is 0.864."),
  "GPW-T-WSR1-POS-B009-O01--r01" = list(method = "t_wilcoxon_signed_rank_one_sample", mode = "post_hoc", branch = "Lehmann d",
    inputs = list(d = 0.1, alpha = 0.05, n_total = 649),
    options = list(method = "lehmann", lehmann_input_method = "predefined_parent", effect_size_parameter = "d", continuity_correction = "off", distribution_plot = "normalized_z", tails = "one", parent_distribution = "normal"),
    outputs = c("Critical z" = "1.6448536", "Power (1-\u03b2 err prob)" = "0.8000589", "Moment p1" = "0.5398278", "Moment p2" = "0.5562315", "Moment p3" = "0.3913924")),
  "SUP-G10-T-WSR1-POS-S01--r01" = list(method = "t_wilcoxon_signed_rank_one_sample", mode = "post_hoc", branch = "Lehmann moments",
    inputs = list(p3 = 0.68, p2 = 0.8, p1 = 0.7, alpha = 0.05, n_total = 20),
    options = list(method = "lehmann", lehmann_input_method = "manual_moments", continuity_correction = "off", distribution_plot = "normalized_z", tails = "one"),
    outputs = c("Critical z" = "1.6448536", "Power (1-\u03b2 err prob)" = "0.8098621")),
  "GPW-T-WSR1-POS-B053-O01--r01" = list(method = "t_wilcoxon_signed_rank_one_sample", mode = "post_hoc", branch = "Lehmann p1 +cc",
    inputs = list(relative_sd = 1, p1 = 0.6382, alpha = 0.05, n_total = 50),
    options = list(method = "lehmann", lehmann_input_method = "predefined_parent", effect_size_parameter = "p1", continuity_correction = "on", distribution_plot = "normalized_z", tails = "one", parent_distribution = "normal"),
    outputs = c("Critical z" = "1.6448536", "Power (1-\u03b2 err prob)" = "0.7740460", "Effect size d" = "0.3536516", "Moment p2" = "0.6915114", "Moment p3" = "0.5463044")),
  "GPW-T-WSR1-POS-B054-O01--r01" = list(method = "t_wilcoxon_signed_rank_one_sample", mode = "post_hoc", branch = "Lehmann p1",
    inputs = list(relative_sd = 1, p1 = 0.5282, alpha = 0.05, n_total = 649),
    options = list(method = "lehmann", lehmann_input_method = "predefined_parent", effect_size_parameter = "p1", continuity_correction = "off", distribution_plot = "normalized_z", tails = "one", parent_distribution = "normal"),
    outputs = c("Critical z" = "1.6448536", "Power (1-\u03b2 err prob)" = "0.5449705", "Effect size d" = "0.0707459", "Moment p2" = "0.5398476", "Moment p3" = "0.3740982")),
  "GPW-T-WSR1-POS-B065-O01--r01" = list(method = "t_wilcoxon_signed_rank_one_sample", mode = "post_hoc", branch = "Lehmann d +cc /T",
    inputs = list(d = 0.5, alpha = 0.05, n_total = 50),
    options = list(method = "lehmann", lehmann_input_method = "predefined_parent", effect_size_parameter = "d", continuity_correction = "on", distribution_plot = "normal_rank_sum", tails = "one", parent_distribution = "normal"),
    outputs = c("Lower critical T" = "808.3933", "Upper critical T" = "808.3933", "Power (1-\u03b2 err prob)" = "0.9678588", "Moment p1" = "0.6914625", "Moment p2" = "0.7602499", "Moment p3" = "0.6337020")),
  "GPW-T-WSR1-POS-B083-O01--r01" = list(method = "t_wilcoxon_signed_rank_one_sample", mode = "post_hoc", branch = "Lehmann moments +cc /T",
    inputs = list(p3 = 0.8206572, p2 = 0.8890997, p1 = 0.8387046, alpha = 0.05, n_total = 50),
    options = list(method = "lehmann", lehmann_input_method = "manual_moments", continuity_correction = "on", distribution_plot = "normal_rank_sum", tails = "one"),
    outputs = c("Lower critical T" = "808.3933", "Upper critical T" = "808.3933", "Power (1-\u03b2 err prob)" = "0.9999999")),
  "GPW-T-WSR1-SEN-B001-O01--r01" = list(method = "t_wilcoxon_signed_rank_one_sample", mode = "sensitivity", branch = "A.R.E. parent",
    inputs = list(alpha = 0.05, power = 0.8, n_total = 50),
    options = list(method = "are", are_specification = "predefined_parent", tails = "one", parent_distribution = "normal"),
    outputs = c("Noncentrality parameter \u03b4" = "2.5232689", "Critical t" = "1.6781097", "Df" = "46.7464829", "Effect size d" = "0.3651681")),
  "GPW-T-WSR1-SEN-B005-D01--r01" = list(method = "t_wilcoxon_signed_rank_one_sample", mode = "sensitivity", branch = "A.R.E. manual",
    inputs = list(are = 0.2, alpha = 0.001, power = 0.99, n_total = 100000),
    options = list(method = "are", are_specification = "manual_are", tails = "one"),
    outputs = c("Noncentrality parameter \u03b4" = "?", "Critical t" = "?", "Df" = "?", "Effect size d" = "?"),
    message = "The minimal A.R.E. of the Wilcoxon test versus the t test is 0.864."),
  "GPW-T-WSR1-SEN-B008-O01--r01" = list(method = "t_wilcoxon_signed_rank_one_sample", mode = "sensitivity", branch = "Lehmann p1 +cc", target = "p1",
    inputs = list(relative_sd = 1, alpha = 0.05, power = 0.8, n_total = 50),
    options = list(method = "lehmann", lehmann_input_method = "predefined_parent", effect_size_parameter = "p1", continuity_correction = "on", distribution_plot = "normalized_z", tails = "one", parent_distribution = "normal"),
    outputs = c("Critical z" = "1.6448536", "Moment p1" = "0.6428804", "Effect size d" = "0.3661686", "Moment p2" = "0.6977153", "Moment p3" = "0.5539423")),
  "GPW-T-WSR1-SEN-B010-O01--r01" = list(method = "t_wilcoxon_signed_rank_one_sample", mode = "sensitivity", branch = "Lehmann d",
    inputs = list(alpha = 0.05, power = 0.8, n_total = 50),
    options = list(method = "lehmann", lehmann_input_method = "predefined_parent", effect_size_parameter = "d", continuity_correction = "off", distribution_plot = "normalized_z", tails = "one", parent_distribution = "normal"),
    outputs = c("Critical z" = "1.6448536", "Effect size d" = "0.3654205", "Moment p1" = "0.6426012", "Moment p2" = "0.6973460", "Moment p3" = "0.5534864")),
  "GPW-T-WSR1-SEN-B011-O01--r01" = list(method = "t_wilcoxon_signed_rank_one_sample", mode = "sensitivity", branch = "Lehmann p1", target = "p1",
    inputs = list(relative_sd = 1, alpha = 0.05, power = 0.8, n_total = 50),
    options = list(method = "lehmann", lehmann_input_method = "predefined_parent", effect_size_parameter = "p1", continuity_correction = "off", distribution_plot = "normalized_z", tails = "one", parent_distribution = "normal"),
    outputs = c("Critical z" = "1.6448536", "Moment p1" = "0.6426012", "Effect size d" = "0.3654205", "Moment p2" = "0.6973460", "Moment p3" = "0.5534864")),
  "GPW-T-WSR1-SEN-B066-O01--r01" = list(method = "t_wilcoxon_signed_rank_one_sample", mode = "sensitivity", branch = "Lehmann d +cc /T",
    inputs = list(alpha = 0.05, power = 0.8, n_total = 50),
    options = list(method = "lehmann", lehmann_input_method = "predefined_parent", effect_size_parameter = "d", continuity_correction = "on", distribution_plot = "normal_rank_sum", tails = "one", parent_distribution = "normal"),
    outputs = c("Lower critical T" = "808.3933", "Upper critical T" = "808.3933", "Effect size d" = "0.3661686", "Moment p1" = "0.6428804", "Moment p2" = "0.6977153", "Moment p3" = "0.5539423")),
  "GPW-UT-DRAW-UTL-U049-UO01--r01" = list(method = "t_wilcoxon_signed_rank_one_sample", mode = "drawer", drawer_mode = "no selector matched (drawer mode 'from means and SD' not found among the live combo items, radios and group boxes); fields used as displayed",
    drawer_inputs = c("Mean H0" = "10", "Mean H1" = "12", "SD \u03c3" = "4"),
    outputs = c("Effect size d" = "0.5")),
  "GPW-T-WSRP-APR-B001-O01--r01" = list(method = "t_wilcoxon_signed_rank_matched_pairs", mode = "a_priori", branch = "A.R.E. parent",
    inputs = list(dz = 0.5, alpha = 0.05, power = 0.8),
    options = list(method = "are", are_specification = "predefined_parent", tails = "one", parent_distribution = "normal"),
    outputs = c("Noncentrality parameter \u03b4" = "2.5854415", "Critical t" = "1.7062592", "Df" = "25.7380304", "Total sample size" = "28", "Actual power" = "0.8083058")),
  "GPW-T-WSRP-APR-B005-D01--r01" = list(method = "t_wilcoxon_signed_rank_matched_pairs", mode = "a_priori", branch = "A.R.E. manual",
    inputs = list(are = 0.2, dz = 0.02, alpha = 0.001, power = 0.99),
    options = list(method = "are", are_specification = "manual_are", tails = "one"),
    outputs = c("Noncentrality parameter \u03b4" = "?", "Critical t" = "?", "Df" = "?", "Total sample size" = "?", "Actual power" = "?"),
    message = "The minimal A.R.E. of the Wilcoxon test versus the t test is 0.864."),
  "GPW-T-WSRP-APR-B009-O01--r01" = list(method = "t_wilcoxon_signed_rank_matched_pairs", mode = "a_priori", branch = "Lehmann d",
    inputs = list(dz = 0.5, alpha = 0.05, power = 0.8),
    options = list(method = "lehmann", lehmann_input_method = "predefined_parent", effect_size_parameter = "d", continuity_correction = "off", distribution_plot = "normalized_z", tails = "one", parent_distribution = "normal"),
    outputs = c("Critical z" = "1.6448536", "Total sample size" = "28", "Actual power" = "0.8078139", "Moment p1" = "0.6914625", "Moment p2" = "0.7602499", "Moment p3" = "0.6337020")),
  "GPW-T-WSRP-APR-B025-O01--r01" = list(method = "t_wilcoxon_signed_rank_matched_pairs", mode = "a_priori", branch = "Lehmann moments",
    inputs = list(p3 = 0.8206572, p2 = 0.8890997, p1 = 0.8387046, alpha = 0.05, power = 0.8),
    options = list(method = "lehmann", lehmann_input_method = "manual_moments", continuity_correction = "off", distribution_plot = "normalized_z", tails = "one"),
    outputs = c("Critical z" = "1.6448536", "Total sample size" = "11", "Actual power" = "0.8308343")),
  "GPW-T-WSRP-APR-B053-O01--r01" = list(method = "t_wilcoxon_signed_rank_matched_pairs", mode = "a_priori", branch = "Lehmann p1 +cc",
    inputs = list(relative_sd = 1, p1 = 0.6382, alpha = 0.05, power = 0.8),
    options = list(method = "lehmann", lehmann_input_method = "predefined_parent", effect_size_parameter = "p1", continuity_correction = "on", distribution_plot = "normalized_z", tails = "one", parent_distribution = "normal"),
    outputs = c("Critical z" = "1.6448536", "Total sample size" = "54", "Actual power" = "0.8037432", "Effect size dz" = "0.3536516", "Moment p2" = "0.6915114", "Moment p3" = "0.5463044")),
  "GPW-T-WSRP-APR-B054-O01--r01" = list(method = "t_wilcoxon_signed_rank_matched_pairs", mode = "a_priori", branch = "Lehmann p1",
    inputs = list(relative_sd = 1, p1 = 0.6382, alpha = 0.05, power = 0.8),
    options = list(method = "lehmann", lehmann_input_method = "predefined_parent", effect_size_parameter = "p1", continuity_correction = "off", distribution_plot = "normalized_z", tails = "one", parent_distribution = "normal"),
    outputs = c("Critical z" = "1.6448536", "Total sample size" = "54", "Actual power" = "0.8050557", "Effect size dz" = "0.3536516", "Moment p2" = "0.6915114", "Moment p3" = "0.5463044")),
  "GPW-T-WSRP-APR-B065-O01--r01" = list(method = "t_wilcoxon_signed_rank_matched_pairs", mode = "a_priori", branch = "Lehmann d +cc /T",
    inputs = list(dz = 0.5, alpha = 0.05, power = 0.8),
    options = list(method = "lehmann", lehmann_input_method = "predefined_parent", effect_size_parameter = "d", continuity_correction = "on", distribution_plot = "normal_rank_sum", tails = "one", parent_distribution = "normal"),
    outputs = c("Lower critical T" = "275.7332", "Upper critical T" = "275.7332", "Total sample size" = "28", "Actual power" = "0.8040183", "Moment p1" = "0.6914625", "Moment p2" = "0.7602499", "Moment p3" = "0.6337020")),
  "GPW-T-WSRP-APR-B083-O01--r01" = list(method = "t_wilcoxon_signed_rank_matched_pairs", mode = "a_priori", branch = "Lehmann moments +cc /T",
    inputs = list(p3 = 0.8206572, p2 = 0.8890997, p1 = 0.8387046, alpha = 0.05, power = 0.8),
    options = list(method = "lehmann", lehmann_input_method = "manual_moments", continuity_correction = "on", distribution_plot = "normal_rank_sum", tails = "one"),
    outputs = c("Lower critical T" = "52.0000337", "Upper critical T" = "52.0000337", "Total sample size" = "11", "Actual power" = "0.8119807")),
  "GPW-T-WSRP-COM-B001-O01--r01" = list(method = "t_wilcoxon_signed_rank_matched_pairs", mode = "compromise", branch = "A.R.E. parent",
    inputs = list(dz = 0.5, error_ratio = 1, n_total = 50),
    options = list(method = "are", are_specification = "predefined_parent", tails = "one", parent_distribution = "normal"),
    outputs = c("Noncentrality parameter \u03b4" = "3.4549415", "Critical t" = "1.7365483", "Df" = "46.7464829", "\u03b1 err prob" = "0.0445268", "\u03b2 err prob" = "0.0445268", "Power (1-\u03b2 err prob)" = "0.9554732")),
  "GPW-T-WSRP-COM-B005-D01--r01" = list(method = "t_wilcoxon_signed_rank_matched_pairs", mode = "compromise", branch = "A.R.E. manual",
    inputs = list(are = 0.2, dz = 0.02, error_ratio = 0.01, n_total = 100000),
    options = list(method = "are", are_specification = "manual_are", tails = "one"),
    outputs = c("Noncentrality parameter \u03b4" = "?", "Critical t" = "?", "Df" = "?", "\u03b1 err prob" = "?", "\u03b2 err prob" = "?", "Power (1-\u03b2 err prob)" = "?"),
    message = "The minimal A.R.E. of the Wilcoxon test versus the t test is 0.864."),
  "GPW-T-WSRP-COM-B009-O01--r01" = list(method = "t_wilcoxon_signed_rank_matched_pairs", mode = "compromise", branch = "Lehmann d",
    inputs = list(dz = 0.5, error_ratio = 1, n_total = 50),
    options = list(method = "lehmann", lehmann_input_method = "predefined_parent", effect_size_parameter = "d", continuity_correction = "off", distribution_plot = "normalized_z", tails = "one", parent_distribution = "normal"),
    outputs = c("Critical z" = "1.7401320", "\u03b1 err prob" = "0.0409179", "\u03b2 err prob" = "0.0409179", "Power (1-\u03b2 err prob)" = "0.9590821", "Moment p1" = "0.6914625", "Moment p2" = "0.7602499", "Moment p3" = "0.6337020")),
  "GPW-T-WSRP-COM-B025-O01--r01" = list(method = "t_wilcoxon_signed_rank_matched_pairs", mode = "compromise", branch = "Lehmann moments",
    inputs = list(p3 = 0.8206572, p2 = 0.8890997, p1 = 0.8387046, error_ratio = 1, n_total = 50),
    options = list(method = "lehmann", lehmann_input_method = "manual_moments", continuity_correction = "off", distribution_plot = "normalized_z", tails = "one"),
    outputs = c("Critical z" = "2.9694897", "\u03b1 err prob" = "0.001491474", "\u03b2 err prob" = "0.001491474", "Power (1-\u03b2 err prob)" = "0.9985085")),
  "GPW-T-WSRP-COM-B053-O01--r01" = list(method = "t_wilcoxon_signed_rank_matched_pairs", mode = "compromise", branch = "Lehmann p1 +cc",
    inputs = list(relative_sd = 1, p1 = 0.6382, error_ratio = 1, n_total = 50),
    options = list(method = "lehmann", lehmann_input_method = "predefined_parent", effect_size_parameter = "p1", continuity_correction = "on", distribution_plot = "normalized_z", tails = "one", parent_distribution = "normal"),
    outputs = c("Critical z" = "1.2205015", "\u03b1 err prob" = "0.1111374", "\u03b2 err prob" = "0.1111374", "Power (1-\u03b2 err prob)" = "0.8888626", "Effect size dz" = "0.3536516", "Moment p2" = "0.6915114", "Moment p3" = "0.5463044")),
  "GPW-T-WSRP-COM-B054-O01--r01" = list(method = "t_wilcoxon_signed_rank_matched_pairs", mode = "compromise", branch = "Lehmann p1",
    inputs = list(relative_sd = 1, p1 = 0.6382, error_ratio = 1, n_total = 50),
    options = list(method = "lehmann", lehmann_input_method = "predefined_parent", effect_size_parameter = "p1", continuity_correction = "off", distribution_plot = "normalized_z", tails = "one", parent_distribution = "normal"),
    outputs = c("Critical z" = "1.2230335", "\u03b1 err prob" = "0.1106585", "\u03b2 err prob" = "0.1106585", "Power (1-\u03b2 err prob)" = "0.8893415", "Effect size dz" = "0.3536516", "Moment p2" = "0.6915114", "Moment p3" = "0.5463044")),
  "GPW-T-WSRP-COM-B065-O01--r01" = list(method = "t_wilcoxon_signed_rank_matched_pairs", mode = "compromise", branch = "Lehmann d +cc /T",
    inputs = list(dz = 0.5, error_ratio = 1, n_total = 50),
    options = list(method = "lehmann", lehmann_input_method = "predefined_parent", effect_size_parameter = "d", continuity_correction = "on", distribution_plot = "normal_rank_sum", tails = "one", parent_distribution = "normal"),
    outputs = c("Lower critical T" = "817.9888", "Upper critical T" = "817.9888", "\u03b1 err prob" = "0.0411510", "\u03b2 err prob" = "0.0411510", "Power (1-\u03b2 err prob)" = "0.9588490", "Moment p1" = "0.6914625", "Moment p2" = "0.7602499", "Moment p3" = "0.6337020")),
  "GPW-T-WSRP-COM-B083-O01--r01" = list(method = "t_wilcoxon_signed_rank_matched_pairs", mode = "compromise", branch = "Lehmann moments +cc /T",
    inputs = list(p3 = 0.8206572, p2 = 0.8890997, p1 = 0.8387046, error_ratio = 1, n_total = 50),
    options = list(method = "lehmann", lehmann_input_method = "manual_moments", continuity_correction = "on", distribution_plot = "normal_rank_sum", tails = "one"),
    outputs = c("Lower critical T" = "945.303", "Upper critical T" = "945.303", "\u03b1 err prob" = "0.001506142", "\u03b2 err prob" = "0.001506142", "Power (1-\u03b2 err prob)" = "0.9984939")),
  "GPW-T-WSRP-CRI-B001-O01--r01" = list(method = "t_wilcoxon_signed_rank_matched_pairs", mode = "criterion", branch = "A.R.E. parent",
    inputs = list(dz = 0.5, power = 0.8, n_total = 50),
    options = list(method = "are", are_specification = "predefined_parent", tails = "one", parent_distribution = "normal"),
    outputs = c("Noncentrality parameter \u03b4" = "3.4549415", "Critical t" = "2.5975827", "Df" = "46.7464829", "\u03b1 err prob" = "0.006257226")),
  "GPW-T-WSRP-CRI-B005-D01--r01" = list(method = "t_wilcoxon_signed_rank_matched_pairs", mode = "criterion", branch = "A.R.E. manual",
    inputs = list(are = 0.2, dz = 0.02, power = 0.9999, n_total = 100000),
    options = list(method = "are", are_specification = "manual_are", tails = "one"),
    outputs = c("Noncentrality parameter \u03b4" = "?", "Critical t" = "?", "Df" = "?", "\u03b1 err prob" = "?"),
    message = "The minimal A.R.E. of the Wilcoxon test versus the t test is 0.864."),
  "GPW-T-WSRP-CRI-B015-O01--r01" = list(method = "t_wilcoxon_signed_rank_matched_pairs", mode = "criterion", branch = "Lehmann d",
    inputs = list(dz = 0.5, power = 0.8, n_total = 50),
    options = list(method = "lehmann", lehmann_input_method = "predefined_parent", effect_size_parameter = "d", continuity_correction = "off", distribution_plot = "normalized_z", tails = "one", parent_distribution = "logistic"),
    outputs = c("Critical z" = "2.6432172", "\u03b1 err prob" = "0.004106117", "Moment p1" = "0.7123653", "Moment p2" = "0.7726573", "Moment p3" = "0.6515369")),
  "GPW-T-WSRP-CRI-B025-O01--r01" = list(method = "t_wilcoxon_signed_rank_matched_pairs", mode = "criterion", branch = "Lehmann moments",
    inputs = list(p3 = 0.8206572, p2 = 0.8890997, p1 = 0.8387046, power = 0.8, n_total = 50),
    options = list(method = "lehmann", lehmann_input_method = "manual_moments", continuity_correction = "off", distribution_plot = "normalized_z", tails = "one"),
    outputs = c("Critical z" = "4.2558879", "\u03b1 err prob" = "1.041105e-05")),
  "GPW-T-WSRP-CRI-B053-O01--r01" = list(method = "t_wilcoxon_signed_rank_matched_pairs", mode = "criterion", branch = "Lehmann p1 +cc",
    inputs = list(relative_sd = 1, p1 = 0.6382, power = 0.8, n_total = 50),
    options = list(method = "lehmann", lehmann_input_method = "predefined_parent", effect_size_parameter = "p1", continuity_correction = "on", distribution_plot = "normalized_z", tails = "one", parent_distribution = "normal"),
    outputs = c("Critical z" = "1.5638522", "\u03b1 err prob" = "0.0589261", "Effect size dz" = "0.3536516", "Moment p2" = "0.6915114", "Moment p3" = "0.5463044")),
  "GPW-T-WSRP-CRI-B054-O01--r01" = list(method = "t_wilcoxon_signed_rank_matched_pairs", mode = "criterion", branch = "Lehmann p1",
    inputs = list(relative_sd = 1, p1 = 0.6382, power = 0.8, n_total = 50),
    options = list(method = "lehmann", lehmann_input_method = "predefined_parent", effect_size_parameter = "p1", continuity_correction = "off", distribution_plot = "normalized_z", tails = "one", parent_distribution = "normal"),
    outputs = c("Critical z" = "1.5686788", "\u03b1 err prob" = "0.0583614", "Effect size dz" = "0.3536516", "Moment p2" = "0.6915114", "Moment p3" = "0.5463044")),
  "GPW-T-WSRP-CRI-B065-O01--r01" = list(method = "t_wilcoxon_signed_rank_matched_pairs", mode = "criterion", branch = "Lehmann d +cc /T",
    inputs = list(dz = 0.5, power = 0.8, n_total = 50),
    options = list(method = "lehmann", lehmann_input_method = "predefined_parent", effect_size_parameter = "d", continuity_correction = "on", distribution_plot = "normal_rank_sum", tails = "one", parent_distribution = "normal"),
    outputs = c("Lower critical T" = "894.2425", "Upper critical T" = "894.2425", "\u03b1 err prob" = "0.006688345", "Moment p1" = "0.6914625", "Moment p2" = "0.7602499", "Moment p3" = "0.6337020")),
  "GPW-T-WSRP-CRI-B083-O01--r01" = list(method = "t_wilcoxon_signed_rank_matched_pairs", mode = "criterion", branch = "Lehmann moments +cc /T",
    inputs = list(p3 = 0.8206572, p2 = 0.8890997, p1 = 0.8387046, power = 0.8, n_total = 50),
    options = list(method = "lehmann", lehmann_input_method = "manual_moments", continuity_correction = "on", distribution_plot = "normal_rank_sum", tails = "one"),
    outputs = c("Lower critical T" = "1078.375", "Upper critical T" = "1078.375", "\u03b1 err prob" = "1.0638e-05")),
  "GPW-T-WSRP-POS-B001-O01--r01" = list(method = "t_wilcoxon_signed_rank_matched_pairs", mode = "post_hoc", branch = "A.R.E. parent",
    inputs = list(dz = 0.5, alpha = 0.05, n_total = 50),
    options = list(method = "are", are_specification = "predefined_parent", tails = "one", parent_distribution = "normal"),
    outputs = c("Noncentrality parameter \u03b4" = "3.4549415", "Critical t" = "1.6781097", "Df" = "46.7464829", "Power (1-\u03b2 err prob)" = "0.9607484")),
  "GPW-T-WSRP-POS-B005-D01--r01" = list(method = "t_wilcoxon_signed_rank_matched_pairs", mode = "post_hoc", branch = "A.R.E. manual",
    inputs = list(are = 0.2, dz = 0.02, alpha = 1e-04, n_total = 100000),
    options = list(method = "are", are_specification = "manual_are", tails = "one"),
    outputs = c("Noncentrality parameter \u03b4" = "?", "Critical t" = "?", "Df" = "?", "Power (1-\u03b2 err prob)" = "?"),
    message = "The minimal A.R.E. of the Wilcoxon test versus the t test is 0.864."),
  "GPW-T-WSRP-POS-B009-O01--r01" = list(method = "t_wilcoxon_signed_rank_matched_pairs", mode = "post_hoc", branch = "Lehmann d",
    inputs = list(dz = 0.5, alpha = 0.05, n_total = 50),
    options = list(method = "lehmann", lehmann_input_method = "predefined_parent", effect_size_parameter = "d", continuity_correction = "off", distribution_plot = "normalized_z", tails = "one", parent_distribution = "normal"),
    outputs = c("Critical z" = "1.6448536", "Power (1-\u03b2 err prob)" = "0.9682796", "Moment p1" = "0.6914625", "Moment p2" = "0.7602499", "Moment p3" = "0.6337020")),
  "GPW-T-WSRP-POS-B025-O01--r01" = list(method = "t_wilcoxon_signed_rank_matched_pairs", mode = "post_hoc", branch = "Lehmann moments",
    inputs = list(p3 = 0.8206572, p2 = 0.8890997, p1 = 0.8387046, alpha = 0.05, n_total = 50),
    options = list(method = "lehmann", lehmann_input_method = "manual_moments", continuity_correction = "off", distribution_plot = "normalized_z", tails = "one"),
    outputs = c("Critical z" = "1.6448536", "Power (1-\u03b2 err prob)" = "0.9999999")),
  "GPW-T-WSRP-POS-B053-O01--r01" = list(method = "t_wilcoxon_signed_rank_matched_pairs", mode = "post_hoc", branch = "Lehmann p1 +cc",
    inputs = list(relative_sd = 1, p1 = 0.6382, alpha = 0.05, n_total = 50),
    options = list(method = "lehmann", lehmann_input_method = "predefined_parent", effect_size_parameter = "p1", continuity_correction = "on", distribution_plot = "normalized_z", tails = "one", parent_distribution = "normal"),
    outputs = c("Critical z" = "1.6448536", "Power (1-\u03b2 err prob)" = "0.7740460", "Effect size dz" = "0.3536516", "Moment p2" = "0.6915114", "Moment p3" = "0.5463044")),
  "GPW-T-WSRP-POS-B054-O01--r01" = list(method = "t_wilcoxon_signed_rank_matched_pairs", mode = "post_hoc", branch = "Lehmann p1",
    inputs = list(relative_sd = 1, p1 = 0.6382, alpha = 0.05, n_total = 50),
    options = list(method = "lehmann", lehmann_input_method = "predefined_parent", effect_size_parameter = "p1", continuity_correction = "off", distribution_plot = "normalized_z", tails = "one", parent_distribution = "normal"),
    outputs = c("Critical z" = "1.6448536", "Power (1-\u03b2 err prob)" = "0.7756440", "Effect size dz" = "0.3536516", "Moment p2" = "0.6915114", "Moment p3" = "0.5463044")),
  "GPW-T-WSRP-POS-B065-O01--r01" = list(method = "t_wilcoxon_signed_rank_matched_pairs", mode = "post_hoc", branch = "Lehmann d +cc /T",
    inputs = list(dz = 0.5, alpha = 0.05, n_total = 50),
    options = list(method = "lehmann", lehmann_input_method = "predefined_parent", effect_size_parameter = "d", continuity_correction = "on", distribution_plot = "normal_rank_sum", tails = "one", parent_distribution = "normal"),
    outputs = c("Lower critical T" = "808.3933", "Upper critical T" = "808.3933", "Power (1-\u03b2 err prob)" = "0.9678588", "Moment p1" = "0.6914625", "Moment p2" = "0.7602499", "Moment p3" = "0.6337020")),
  "GPW-T-WSRP-POS-B083-O01--r01" = list(method = "t_wilcoxon_signed_rank_matched_pairs", mode = "post_hoc", branch = "Lehmann moments +cc /T",
    inputs = list(p3 = 0.8206572, p2 = 0.8890997, p1 = 0.8387046, alpha = 0.05, n_total = 50),
    options = list(method = "lehmann", lehmann_input_method = "manual_moments", continuity_correction = "on", distribution_plot = "normal_rank_sum", tails = "one"),
    outputs = c("Lower critical T" = "808.3933", "Upper critical T" = "808.3933", "Power (1-\u03b2 err prob)" = "0.9999999")),
  "GPW-T-WSRP-SEN-B001-O01--r01" = list(method = "t_wilcoxon_signed_rank_matched_pairs", mode = "sensitivity", branch = "A.R.E. parent",
    inputs = list(alpha = 0.05, power = 0.8, n_total = 50),
    options = list(method = "are", are_specification = "predefined_parent", tails = "one", parent_distribution = "normal"),
    outputs = c("Noncentrality parameter \u03b4" = "2.5232689", "Critical t" = "1.6781097", "Df" = "46.7464829", "Effect size dz" = "0.3651681")),
  "GPW-T-WSRP-SEN-B005-D01--r01" = list(method = "t_wilcoxon_signed_rank_matched_pairs", mode = "sensitivity", branch = "A.R.E. manual",
    inputs = list(are = 0.2, alpha = 0.001, power = 0.99, n_total = 100000),
    options = list(method = "are", are_specification = "manual_are", tails = "one"),
    outputs = c("Noncentrality parameter \u03b4" = "?", "Critical t" = "?", "Df" = "?", "Effect size dz" = "?"),
    message = "The minimal A.R.E. of the Wilcoxon test versus the t test is 0.864."),
  "GPW-T-WSRP-SEN-B008-O01--r01" = list(method = "t_wilcoxon_signed_rank_matched_pairs", mode = "sensitivity", branch = "Lehmann p1 +cc", target = "p1",
    inputs = list(relative_sd = 1, alpha = 0.05, power = 0.8, n_total = 50),
    options = list(method = "lehmann", lehmann_input_method = "predefined_parent", effect_size_parameter = "p1", continuity_correction = "on", distribution_plot = "normalized_z", tails = "one", parent_distribution = "normal"),
    outputs = c("Critical z" = "1.6448536", "Moment p1" = "0.6428804", "Effect size dz" = "0.3661686", "Moment p2" = "0.6977153", "Moment p3" = "0.5539423")),
  "GPW-T-WSRP-SEN-B010-O01--r01" = list(method = "t_wilcoxon_signed_rank_matched_pairs", mode = "sensitivity", branch = "Lehmann d",
    inputs = list(alpha = 0.05, power = 0.8, n_total = 50),
    options = list(method = "lehmann", lehmann_input_method = "predefined_parent", effect_size_parameter = "d", continuity_correction = "off", distribution_plot = "normalized_z", tails = "one", parent_distribution = "normal"),
    outputs = c("Critical z" = "1.6448536", "Effect size dz" = "0.3654205", "Moment p1" = "0.6426012", "Moment p2" = "0.6973460", "Moment p3" = "0.5534864")),
  "GPW-T-WSRP-SEN-B011-O01--r01" = list(method = "t_wilcoxon_signed_rank_matched_pairs", mode = "sensitivity", branch = "Lehmann p1", target = "p1",
    inputs = list(relative_sd = 1, alpha = 0.05, power = 0.8, n_total = 50),
    options = list(method = "lehmann", lehmann_input_method = "predefined_parent", effect_size_parameter = "p1", continuity_correction = "off", distribution_plot = "normalized_z", tails = "one", parent_distribution = "normal"),
    outputs = c("Critical z" = "1.6448536", "Moment p1" = "0.6426012", "Effect size dz" = "0.3654205", "Moment p2" = "0.6973460", "Moment p3" = "0.5534864")),
  "GPW-T-WSRP-SEN-B066-O01--r01" = list(method = "t_wilcoxon_signed_rank_matched_pairs", mode = "sensitivity", branch = "Lehmann d +cc /T",
    inputs = list(alpha = 0.05, power = 0.8, n_total = 50),
    options = list(method = "lehmann", lehmann_input_method = "predefined_parent", effect_size_parameter = "d", continuity_correction = "on", distribution_plot = "normal_rank_sum", tails = "one", parent_distribution = "normal"),
    outputs = c("Lower critical T" = "808.3933", "Upper critical T" = "808.3933", "Effect size dz" = "0.3661686", "Moment p1" = "0.6428804", "Moment p2" = "0.6977153", "Moment p3" = "0.5539423")),
  "GPW-UT-DRAW-UTL-U050-UO01--r01" = list(method = "t_wilcoxon_signed_rank_matched_pairs", mode = "drawer", drawer_mode = "group 'From differences'",
    drawer_inputs = c("Mean of difference" = "0.5", "SD of difference" = "1"),
    outputs = c("Effect size dz" = "0.5")),
  "GPW-UT-DRAW-UTL-U051-UO01--r01" = list(method = "t_wilcoxon_signed_rank_matched_pairs", mode = "drawer", drawer_mode = "group 'From group parameters'",
    drawer_inputs = c("Mean group 1" = "6500", "Mean group 2" = "5600", "SD group 1" = "1500", "SD group 2" = "1300", "Correlation between groups" = "0.85"),
    outputs = c("Effect size dz" = "1.13842")),
  "GPW-T-WMW-APR-B001-O01--r01" = list(method = "t_wilcoxon_mann_whitney", mode = "a_priori", branch = "A.R.E. parent",
    inputs = list(d = 0.5, alpha = 0.05, power = 0.8, allocation_ratio = 1),
    options = list(method = "are", are_specification = "predefined_parent", tails = "one", parent_distribution = "normal"),
    outputs = c("Noncentrality parameter \u03b4" = "2.5152354", "Critical t" = "1.6603560", "Df" = "99.2225438", "Sample size group 1" = "53", "Sample size group 2" = "53", "Total sample size" = "106", "Actual power" = "0.8032180")),
  "GPW-T-WMW-APR-B005-O01--r01" = list(method = "t_wilcoxon_mann_whitney", mode = "a_priori", branch = "A.R.E. manual",
    inputs = list(are = 0.955, d = 0.5, alpha = 0.05, power = 0.8, allocation_ratio = 1),
    options = list(method = "are", are_specification = "manual_are", tails = "one"),
    outputs = c("Noncentrality parameter \u03b4" = "2.5153280", "Critical t" = "1.6603548", "Df" = "99.2300000", "Sample size group 1" = "53", "Sample size group 2" = "53", "Total sample size" = "106", "Actual power" = "0.8032439")),
  "GPW-T-WMW-APR-B015-O01--r01" = list(method = "t_wilcoxon_mann_whitney", mode = "a_priori", branch = "Lehmann d",
    inputs = list(d = 0.5, alpha = 0.05, power = 0.8, allocation_ratio = 1),
    options = list(method = "lehmann", lehmann_input_method = "predefined_parent", effect_size_parameter = "d", continuity_correction = "off", distribution_plot = "normalized_z", tails = "one", parent_distribution = "logistic"),
    outputs = c("Critical z" = "1.6448536", "Sample size group 1" = "47", "Sample size group 2" = "47", "Total sample size" = "94", "Actual power" = "0.8045264", "Moment p1" = "0.6471246", "Moment p2" = "0.4935195")),
  "GPW-T-WMW-APR-B054-O01--r01" = list(method = "t_wilcoxon_mann_whitney", mode = "a_priori", branch = "Lehmann p1",
    inputs = list(relative_sd = 1, p1 = 0.6382, alpha = 0.05, power = 0.8, allocation_ratio = 1),
    options = list(method = "lehmann", lehmann_input_method = "predefined_parent", effect_size_parameter = "p1", continuity_correction = "off", distribution_plot = "normalized_z", tails = "one", parent_distribution = "normal"),
    outputs = c("Critical z" = "1.6448536", "Sample size group 1" = "53", "Sample size group 2" = "53", "Total sample size" = "106", "Actual power" = "0.8022293", "Effect size d" = "0.5001389", "Moment p2" = "0.4826356")),
  "GPW-T-WMW-APR-B055-O01--r01" = list(method = "t_wilcoxon_mann_whitney", mode = "a_priori", branch = "Lehmann p1 +cc",
    inputs = list(relative_sd = 1, p1 = 0.6382, alpha = 0.05, power = 0.8, allocation_ratio = 1),
    options = list(method = "lehmann", lehmann_input_method = "predefined_parent", effect_size_parameter = "p1", continuity_correction = "on", distribution_plot = "normalized_z", tails = "one", parent_distribution = "logistic"),
    outputs = c("Critical z" = "1.6448536", "Sample size group 1" = "53", "Sample size group 2" = "53", "Total sample size" = "106", "Actual power" = "0.8006750", "Effect size d" = "0.4681291", "Moment p2" = "0.4830505")),
  "GPW-T-WMW-APR-B065-O01--r01" = list(method = "t_wilcoxon_mann_whitney", mode = "a_priori", branch = "Lehmann d +cc /T",
    inputs = list(d = 0.5, alpha = 0.05, power = 0.8, allocation_ratio = 1),
    options = list(method = "lehmann", lehmann_input_method = "predefined_parent", effect_size_parameter = "d", continuity_correction = "on", distribution_plot = "normal_rank_sum", tails = "one", parent_distribution = "normal"),
    outputs = c("Lower critical T" = "1665.318", "Upper critical T" = "1665.318", "Sample size group 1" = "53", "Sample size group 2" = "53", "Total sample size" = "106", "Actual power" = "0.8011061", "Moment p1" = "0.6381632", "Moment p2" = "0.4825929")),
  "GPW-T-WMW-COM-B001-O01--r01" = list(method = "t_wilcoxon_mann_whitney", mode = "compromise", branch = "A.R.E. parent",
    inputs = list(d = 0.5, error_ratio = 1, n1 = 50, n2 = 50),
    options = list(method = "are", are_specification = "predefined_parent", tails = "one", parent_distribution = "normal"),
    outputs = c("Noncentrality parameter \u03b4" = "2.4430126", "Critical t" = "1.2247727", "Df" = "93.4929659", "\u03b1 err prob" = "0.1118689", "\u03b2 err prob" = "0.1118689", "Power (1-\u03b2 err prob)" = "0.8881311")),
  "GPW-T-WMW-COM-B005-O01--r01" = list(method = "t_wilcoxon_mann_whitney", mode = "compromise", branch = "A.R.E. manual",
    inputs = list(are = 0.955, d = 0.5, error_ratio = 1, n1 = 50, n2 = 50),
    options = list(method = "are", are_specification = "manual_are", tails = "one"),
    outputs = c("Noncentrality parameter \u03b4" = "2.4431025", "Critical t" = "1.2248176", "Df" = "93.5000000", "\u03b1 err prob" = "0.1118604", "\u03b2 err prob" = "0.1118604", "Power (1-\u03b2 err prob)" = "0.8881396")),
  "GPW-T-WMW-COM-B015-O01--r01" = list(method = "t_wilcoxon_mann_whitney", mode = "compromise", branch = "Lehmann d",
    inputs = list(d = 0.5, error_ratio = 1, n1 = 50, n2 = 50),
    options = list(method = "lehmann", lehmann_input_method = "predefined_parent", effect_size_parameter = "d", continuity_correction = "off", distribution_plot = "normalized_z", tails = "one", parent_distribution = "logistic"),
    outputs = c("Critical z" = "1.3020875", "\u03b1 err prob" = "0.0964432", "\u03b2 err prob" = "0.0964432", "Power (1-\u03b2 err prob)" = "0.9035568", "Moment p1" = "0.6471246", "Moment p2" = "0.4935195")),
  "GPW-T-WMW-COM-B055-O01--r01" = list(method = "t_wilcoxon_mann_whitney", mode = "compromise", branch = "Lehmann p1 +cc",
    inputs = list(relative_sd = 1, p1 = 0.6382, error_ratio = 1, n1 = 50, n2 = 50),
    options = list(method = "lehmann", lehmann_input_method = "predefined_parent", effect_size_parameter = "p1", continuity_correction = "on", distribution_plot = "normalized_z", tails = "one", parent_distribution = "logistic"),
    outputs = c("Critical z" = "1.2173971", "\u03b1 err prob" = "0.1117266", "\u03b2 err prob" = "0.1117266", "Power (1-\u03b2 err prob)" = "0.8882734", "Effect size d" = "0.4681291", "Moment p2" = "0.4830505")),
  "GPW-T-WMW-COM-B054-O01--r01" = list(method = "t_wilcoxon_mann_whitney", mode = "compromise", branch = "Lehmann p1",
    inputs = list(relative_sd = 1, p1 = 0.6382, error_ratio = 1, n1 = 50, n2 = 50),
    options = list(method = "lehmann", lehmann_input_method = "predefined_parent", effect_size_parameter = "p1", continuity_correction = "off", distribution_plot = "normalized_z", tails = "one", parent_distribution = "normal"),
    outputs = c("Critical z" = "1.2207461", "\u03b1 err prob" = "0.1110911", "\u03b2 err prob" = "0.1110911", "Power (1-\u03b2 err prob)" = "0.8889089", "Effect size d" = "0.5001389", "Moment p2" = "0.4826356")),
  "GPW-T-WMW-COM-B065-O01--r01" = list(method = "t_wilcoxon_mann_whitney", mode = "compromise", branch = "Lehmann d +cc /T",
    inputs = list(d = 0.5, error_ratio = 1, n1 = 50, n2 = 50),
    options = list(method = "lehmann", lehmann_input_method = "predefined_parent", effect_size_parameter = "d", continuity_correction = "on", distribution_plot = "normal_rank_sum", tails = "one", parent_distribution = "normal"),
    outputs = c("Lower critical T" = "1427.273", "Upper critical T" = "1427.273", "\u03b1 err prob" = "0.1114908", "\u03b2 err prob" = "0.1114908", "Power (1-\u03b2 err prob)" = "0.8885092", "Moment p1" = "0.6381632", "Moment p2" = "0.4825929")),
  "GPW-T-WMW-CRI-B001-O01--r01" = list(method = "t_wilcoxon_mann_whitney", mode = "criterion", branch = "A.R.E. parent",
    inputs = list(d = 0.5, power = 0.8, n1 = 50, n2 = 50),
    options = list(method = "are", are_specification = "predefined_parent", tails = "one", parent_distribution = "normal"),
    outputs = c("Noncentrality parameter \u03b4" = "2.4430126", "Critical t" = "1.5999435", "Df" = "93.4929659", "\u03b1 err prob" = "0.0564925")),
  "GPW-T-WMW-CRI-B005-O01--r01" = list(method = "t_wilcoxon_mann_whitney", mode = "criterion", branch = "A.R.E. manual",
    inputs = list(are = 0.955, d = 0.5, power = 0.8, n1 = 50, n2 = 50),
    options = list(method = "are", are_specification = "manual_are", tails = "one"),
    outputs = c("Noncentrality parameter \u03b4" = "2.4431025", "Critical t" = "1.6000332", "Df" = "93.5000000", "\u03b1 err prob" = "0.0564825")),
  "GPW-T-WMW-CRI-B015-O01--r01" = list(method = "t_wilcoxon_mann_whitney", mode = "criterion", branch = "Lehmann d",
    inputs = list(d = 0.5, power = 0.8, n1 = 50, n2 = 50),
    options = list(method = "lehmann", lehmann_input_method = "predefined_parent", effect_size_parameter = "d", continuity_correction = "off", distribution_plot = "normalized_z", tails = "one", parent_distribution = "logistic"),
    outputs = c("Critical z" = "1.7383121", "\u03b1 err prob" = "0.0410779", "Moment p1" = "0.6471246", "Moment p2" = "0.4935195")),
  "GPW-T-WMW-CRI-B055-O01--r01" = list(method = "t_wilcoxon_mann_whitney", mode = "criterion", branch = "Lehmann p1 +cc",
    inputs = list(relative_sd = 1, p1 = 0.6382, power = 0.8, n1 = 50, n2 = 50),
    options = list(method = "lehmann", lehmann_input_method = "predefined_parent", effect_size_parameter = "p1", continuity_correction = "on", distribution_plot = "normalized_z", tails = "one", parent_distribution = "logistic"),
    outputs = c("Critical z" = "1.5757558", "\u03b1 err prob" = "0.0575411", "Effect size d" = "0.4681291", "Moment p2" = "0.4830505")),
  "GPW-T-WMW-CRI-B056-O01--r01" = list(method = "t_wilcoxon_mann_whitney", mode = "criterion", branch = "Lehmann p1",
    inputs = list(relative_sd = 1, p1 = 0.6382, power = 0.8, n1 = 50, n2 = 50),
    options = list(method = "lehmann", lehmann_input_method = "predefined_parent", effect_size_parameter = "p1", continuity_correction = "off", distribution_plot = "normalized_z", tails = "one", parent_distribution = "logistic"),
    outputs = c("Critical z" = "1.5792027", "\u03b1 err prob" = "0.0571448", "Effect size d" = "0.4681291", "Moment p2" = "0.4830505")),
  "GPW-T-WMW-CRI-B065-O01--r01" = list(method = "t_wilcoxon_mann_whitney", mode = "criterion", branch = "Lehmann d +cc /T",
    inputs = list(d = 0.5, power = 0.8, n1 = 50, n2 = 50),
    options = list(method = "lehmann", lehmann_input_method = "predefined_parent", effect_size_parameter = "d", continuity_correction = "on", distribution_plot = "normal_rank_sum", tails = "one", parent_distribution = "normal"),
    outputs = c("Lower critical T" = "1479.29", "Upper critical T" = "1479.29", "\u03b1 err prob" = "0.0573709", "Moment p1" = "0.6381632", "Moment p2" = "0.4825929")),
  "SUP-G10-T-WMW-POS-D03--r01" = list(method = "t_wilcoxon_mann_whitney", mode = "post_hoc", branch = "A.R.E. parent",
    inputs = list(d = 0.5, alpha = 0.05, n1 = 1, n2 = 2),
    options = list(method = "are", are_specification = "predefined_parent", tails = "one", parent_distribution = "normal"),
    outputs = c("Noncentrality parameter \u03b4" = "0.3989423", "Critical t" = "8.2783370", "Df" = "0.8647890", "Power (1-\u03b2 err prob)" = "0.0771681")),
  "GPW-T-WMW-POS-B005-O01--r01" = list(method = "t_wilcoxon_mann_whitney", mode = "post_hoc", branch = "A.R.E. manual",
    inputs = list(are = 0.955, d = 0.5, alpha = 0.05, n1 = 50, n2 = 50),
    options = list(method = "are", are_specification = "manual_are", tails = "one"),
    outputs = c("Noncentrality parameter \u03b4" = "2.4431025", "Critical t" = "1.6613143", "Df" = "93.5000000", "Power (1-\u03b2 err prob)" = "0.7824544")),
  "GPW-T-WMW-POS-B009-O01--r01" = list(method = "t_wilcoxon_mann_whitney", mode = "post_hoc", branch = "Lehmann d",
    inputs = list(d = 0.5, alpha = 0.05, n1 = 50, n2 = 50),
    options = list(method = "lehmann", lehmann_input_method = "predefined_parent", effect_size_parameter = "d", continuity_correction = "off", distribution_plot = "normalized_z", tails = "one", parent_distribution = "normal"),
    outputs = c("Critical z" = "1.6448536", "Power (1-\u03b2 err prob)" = "0.7805794", "Moment p1" = "0.6381632", "Moment p2" = "0.4825929")),
  "GPW-T-WMW-POS-B053-O01--r01" = list(method = "t_wilcoxon_mann_whitney", mode = "post_hoc", branch = "Lehmann p1 +cc",
    inputs = list(relative_sd = 1, p1 = 0.6382, alpha = 0.05, n1 = 50, n2 = 50),
    options = list(method = "lehmann", lehmann_input_method = "predefined_parent", effect_size_parameter = "p1", continuity_correction = "on", distribution_plot = "normalized_z", tails = "one", parent_distribution = "normal"),
    outputs = c("Critical z" = "1.6448536", "Power (1-\u03b2 err prob)" = "0.7797105", "Effect size d" = "0.5001389", "Moment p2" = "0.4826356")),
  "SUP-G10-T-WMW-POS-P01--r01" = list(method = "t_wilcoxon_mann_whitney", mode = "post_hoc", branch = "Lehmann p1",
    inputs = list(relative_sd = 1, p1 = 0.6962311, alpha = 0.05, n1 = 50, n2 = 50),
    options = list(method = "lehmann", lehmann_input_method = "predefined_parent", effect_size_parameter = "p1", continuity_correction = "off", distribution_plot = "normalized_z", tails = "one", parent_distribution = "normal"),
    outputs = c("Critical z" = "1.6448536", "Power (1-\u03b2 err prob)" = "0.9732339", "Effect size d" = "0.7263277", "Moment p2" = "0.5521107")),
  "GPW-T-WMW-POS-B065-O01--r01" = list(method = "t_wilcoxon_mann_whitney", mode = "post_hoc", branch = "Lehmann d +cc /T",
    inputs = list(d = 0.5, alpha = 0.05, n1 = 50, n2 = 50),
    options = list(method = "lehmann", lehmann_input_method = "predefined_parent", effect_size_parameter = "d", continuity_correction = "on", distribution_plot = "normal_rank_sum", tails = "one", parent_distribution = "normal"),
    outputs = c("Lower critical T" = "1489.098", "Upper critical T" = "1489.098", "Power (1-\u03b2 err prob)" = "0.7795065", "Moment p1" = "0.6381632", "Moment p2" = "0.4825929")),
  "GPW-T-WMW-SEN-B001-O01--r01" = list(method = "t_wilcoxon_mann_whitney", mode = "sensitivity", branch = "A.R.E. parent",
    inputs = list(alpha = 0.05, power = 0.8, n1 = 50, n2 = 50),
    options = list(method = "are", are_specification = "predefined_parent", tails = "one", parent_distribution = "normal"),
    outputs = c("Noncentrality parameter \u03b4" = "2.5046661", "Critical t" = "1.6613155", "Df" = "93.4929659", "Effect size d" = "0.5126183")),
  "GPW-T-WMW-SEN-B005-O01--r01" = list(method = "t_wilcoxon_mann_whitney", mode = "sensitivity", branch = "A.R.E. manual",
    inputs = list(are = 0.955, alpha = 0.05, power = 0.8, n1 = 50, n2 = 50),
    options = list(method = "are", are_specification = "manual_are", tails = "one"),
    outputs = c("Noncentrality parameter \u03b4" = "2.5046647", "Critical t" = "1.6613143", "Df" = "93.5000000", "Effect size d" = "0.5125992")),
  "GPW-T-WMW-SEN-B016-O01--r01" = list(method = "t_wilcoxon_mann_whitney", mode = "sensitivity", branch = "Lehmann d",
    inputs = list(alpha = 0.05, power = 0.8, n1 = 50, n2 = 50),
    options = list(method = "lehmann", lehmann_input_method = "predefined_parent", effect_size_parameter = "d", continuity_correction = "off", distribution_plot = "normalized_z", tails = "one", parent_distribution = "logistic"),
    outputs = c("Critical z" = "1.6448536", "Effect size d" = "0.4812496", "Moment p1" = "0.6418848", "Moment p2" = "0.4873610")),
  "GPW-T-WMW-SEN-B072-O01--r01" = list(method = "t_wilcoxon_mann_whitney", mode = "sensitivity", branch = "Lehmann d +cc /T",
    inputs = list(alpha = 0.05, power = 0.8, n1 = 50, n2 = 50),
    options = list(method = "lehmann", lehmann_input_method = "predefined_parent", effect_size_parameter = "d", continuity_correction = "on", distribution_plot = "normal_rank_sum", tails = "one", parent_distribution = "logistic"),
    outputs = c("Lower critical T" = "1489.098", "Upper critical T" = "1489.098", "Effect size d" = "0.4819396", "Moment p1" = "0.6420782", "Moment p2" = "0.4875877")),
  "GPW-UT-DRAW-UTL-U052-UO01--r02" = list(method = "t_wilcoxon_mann_whitney", mode = "drawer", drawer_mode = "group 'n1 != n2'",
    drawer_inputs = c("Mean group 1" = "0.3", "Mean group 2" = "0.27", "SD \u03c3 within each group" = "0.08"),
    outputs = c("Effect size d" = "0.375")),
  "GPW-UT-DRAW-UTL-U053-UO01--r01" = list(method = "t_wilcoxon_mann_whitney", mode = "drawer", drawer_mode = "group 'n1 = n2'",
    drawer_inputs = c("Mean group 1" = "0", "Mean group 2" = "1", "SD \u03c3 group 1" = "1", "SD \u03c3 group 2" = "2"),
    outputs = c("Effect size d" = "0.6324555")),
  "SUP-G10-T-WMW-APR-M04--r01" = list(method = "t_wilcoxon_mann_whitney", mode = "a_priori", branch = "Lehmann moments",
    inputs = list(p3 = 0.4, p2 = 0.45, p1 = 0.6, alpha = 0.05, power = 0.8, allocation_ratio = 3),
    options = list(method = "lehmann", lehmann_input_method = "manual_moments", continuity_correction = "off", distribution_plot = "normalized_z", tails = "one"),
    outputs = c("Critical z" = "1.6448536", "Sample size group 1" = "68", "Sample size group 2" = "202", "Total sample size" = "270", "Actual power" = "0.8027793")),
  "SUP-G10-T-WMW-COM-M06--r01" = list(method = "t_wilcoxon_mann_whitney", mode = "compromise", branch = "Lehmann moments +cc",
    inputs = list(p3 = 0.44, p2 = 0.4, p1 = 0.58, error_ratio = 2, n1 = 20, n2 = 50),
    options = list(method = "lehmann", lehmann_input_method = "manual_moments", continuity_correction = "on", distribution_plot = "normalized_z", tails = "one"),
    outputs = c("Critical z" = "0.8212607", "\u03b1 err prob" = "0.2057489", "\u03b2 err prob" = "0.4114978", "Power (1-\u03b2 err prob)" = "0.5885022")),
  "SUP-G10-T-WMW-CRI-M05--r01" = list(method = "t_wilcoxon_mann_whitney", mode = "criterion", branch = "Lehmann moments /T",
    inputs = list(p3 = 0.4, p2 = 0.45, p1 = 0.6, power = 0.9, n1 = 30, n2 = 10),
    options = list(method = "lehmann", lehmann_input_method = "manual_moments", continuity_correction = "off", distribution_plot = "normal_rank_sum", tails = "two"),
    outputs = c("Lower critical T" = "143.6604", "Upper critical T" = "156.3396", "\u03b1 err prob" = "0.8430318")),
  "SUP-G10-T-WMW-POS-M01--r01" = list(method = "t_wilcoxon_mann_whitney", mode = "post_hoc", branch = "Lehmann moments",
    inputs = list(p3 = 0.42, p2 = 0.47, p1 = 0.62, alpha = 0.05, n1 = 5, n2 = 40),
    options = list(method = "lehmann", lehmann_input_method = "manual_moments", continuity_correction = "off", distribution_plot = "normalized_z", tails = "one"),
    outputs = c("Critical z" = "1.6448536", "Power (1-\u03b2 err prob)" = "0.2146992")),
  "SUP-G10-T-WMW-POS-M03--r01" = list(method = "t_wilcoxon_mann_whitney", mode = "post_hoc", branch = "Lehmann moments +cc",
    inputs = list(p3 = 0.47, p2 = 0.42, p1 = 0.62, alpha = 0.05, n1 = 12, n2 = 60),
    options = list(method = "lehmann", lehmann_input_method = "manual_moments", continuity_correction = "on", distribution_plot = "normalized_z", tails = "two"),
    outputs = c("Critical z" = "1.9599640", "Power (1-\u03b2 err prob)" = "0.1835562"))
)

g10_exact_labels <- c("Total sample size", "Sample size group 1", "Sample size group 2")

g10_display_ok <- function(ours, shown) {
  if (!is.numeric(ours) || length(ours) != 1L || !is.finite(ours)) return(FALSE)
  m <- regmatches(shown, regexec("^(-?[0-9]*)(\\.([0-9]*))?([eE]([+-]?[0-9]+))?$", shown))[[1]]
  if (!length(m)) return(FALSE)
  dec <- nchar(m[4])
  ex <- if (nzchar(m[6])) as.numeric(m[6]) else 0
  x <- as.numeric(shown)
  half <- 0.5 * 10^(-dec) * 10^ex
  ours >= x - half - 8 * .Machine$double.eps * max(1, abs(x)) && ours <= x + half + 8 * .Machine$double.eps * max(1, abs(x))
}

g10_outputs <- function(res, eff) {
  q <- res$quantities
  s <- res$solved
  a <- res$achieved
  mo <- res$method_outputs
  gs <- res$group_sizes
  lehmann <- identical(mo$power_method, "lehmann")
  crit_t <- if (lehmann) NULL else if (length(q$critical) == 1L) q$critical else
    if (q$ncp < 0) min(q$critical) else max(q$critical)
  d_out <- s[[eff]] %||% mo$effect_size_d
  list("Noncentrality parameter \u03b4" = if (!lehmann) q$ncp, "Critical t" = crit_t, "Df" = if (!lehmann) q$df,
       "Critical z" = if (lehmann) max(mo$critical_z),
       "Lower critical T" = if (lehmann) mo$critical_rank_sum_lower, "Upper critical T" = if (lehmann) mo$critical_rank_sum_upper,
       "Sample size group 1" = if (!is.null(gs)) unname(gs[["n1"]]), "Sample size group 2" = if (!is.null(gs)) unname(gs[["n2"]]),
       "Total sample size" = s$n_total %||% q$n_total, "Actual power" = a$power, "Power (1-\u03b2 err prob)" = a$power,
       "\u03b1 err prob" = s$alpha, "\u03b2 err prob" = a$beta, "Effect size d" = d_out, "Effect size dz" = d_out,
       "Moment p1" = if (lehmann) s$p1 %||% mo$moment_p1, "Moment p2" = if (lehmann) mo$moment_p2,
       "Moment p3" = if (lehmann) mo$moment_p3)
}

g10_drawer_value <- function(cs) {
  num <- function(lab) as.numeric(cs$drawer_inputs[[lab]])
  if (cs$method == "t_wilcoxon_signed_rank_one_sample") {
    return(effect_size(cs$method, "means_sd", mean_h0 = num("Mean H0"), mean_h1 = num("Mean H1"), sd = num("SD \u03c3"))$value)
  }
  if (cs$method == "t_wilcoxon_signed_rank_matched_pairs") {
    if (grepl("From differences", cs$drawer_mode, fixed = TRUE)) {
      return(effect_size(cs$method, "differences", mean_diff = num("Mean of difference"), sd_diff = num("SD of difference"))$value)
    }
    return(effect_size(cs$method, "group_parameters", mean_1 = num("Mean group 1"), mean_2 = num("Mean group 2"),
                       sd_1 = num("SD group 1"), sd_2 = num("SD group 2"), r = num("Correlation between groups"))$magnitude)
  }
  if (grepl("n1 != n2", cs$drawer_mode, fixed = TRUE)) {
    return(effect_size(cs$method, "common_sd", mean_1 = num("Mean group 1"), mean_2 = num("Mean group 2"),
                       sd = num("SD \u03c3 within each group"))$magnitude)
  }
  effect_size(cs$method, "two_sds", mean_1 = num("Mean group 1"), mean_2 = num("Mean group 2"),
              sd_1 = num("SD \u03c3 group 1"), sd_2 = num("SD \u03c3 group 2"))$magnitude
}

g10_check_case <- function(id) {
  cs <- g10_fixture_cases[[id]]
  if (identical(cs$mode, "drawer")) {
    v <- g10_drawer_value(cs)
    for (lab in names(cs$outputs)) expect_true(g10_display_ok(v, cs$outputs[[lab]]), info = paste(id, lab, v))
    return(invisible())
  }
  spec <- tryCatch(analysis_spec(cs$method, cs$mode, inputs = cs$inputs, options = cs$options, target = cs$target),
                   inferguidR_error = function(e) e)
  if (!is.null(cs$message)) {
    # TD-0010: a refusal is matched by a typed refusal without numbers.
    res <- if (inherits(spec, "inferguidR_error")) spec else calculate_power(spec)
    code <- res$status$code
    expect_true(code %in% c("invalid_input", "infeasible"), info = paste(id, code))
    if (!inherits(res, "inferguidR_error")) {
      expect_identical(.ig_finite_paths(unclass(res)[c("solved", "achieved", "quantities")], ""), character(), info = id)
    }
    return(invisible())
  }
  expect_false(inherits(spec, "inferguidR_error"), info = paste(id, if (inherits(spec, "inferguidR_error")) conditionMessage(spec)))
  res <- calculate_power(spec)
  expect_identical(res$status$code, "ok", info = paste(id, res$status$reason))
  outs <- g10_outputs(res, .ig_wilcoxon_info(cs$method)$effect)
  for (lab in names(cs$outputs)) {
    v <- outs[[lab]]
    shown <- cs$outputs[[lab]]
    if (lab %in% g10_exact_labels) {
      expect_identical(as.double(v), as.numeric(shown), info = paste(id, lab))
    } else {
      expect_true(g10_display_ok(v, shown),
                  info = sprintf("%s %s: G*Power %s, engine %.12g", id, lab, shown, if (is.null(v)) NA_real_ else v))
    }
  }
}

for (m in unique(vapply(g10_fixture_cases, `[[`, "", "method"))) {
  local({
    method <- m
    test_that(sprintf("curated G*Power fixtures agree: %s", method), {
      ids <- names(g10_fixture_cases)[vapply(g10_fixture_cases, `[[`, "", "method") == method]
      for (id in ids) g10_check_case(id)
      modes <- unique(vapply(g10_fixture_cases[ids], `[[`, "", "mode"))
      expect_true(all(.ig_get_method(method)$modes %in% modes), info = paste(method, "every mode has a curated case"))
    })
  })
}

test_that("the embedded displayed strings equal the fixture records (development checkout)", {
  root <- normalizePath(testthat::test_path("..", ".."), mustWork = FALSE)
  fx <- file.path(root, "validation", "fixtures")
  skip_if_not(dir.exists(fx), "fixture records are not part of the installed package")
  for (id in names(g10_fixture_cases)) {
    cs <- g10_fixture_cases[[id]]
    f <- c(file.path(fx, "emulator_capture", paste0(id, ".json")),
           file.path(fx, "emulator_capture_supplement", "g10", paste0(id, ".json")))
    f <- f[file.exists(f)]
    expect_length(f, 1L)
    rec <- jsonlite::parse_json(readChar(f, file.size(f), useBytes = TRUE), simplifyVector = FALSE)
    outs <- c(rec$outputs, rec$effect_size_drawer$outputs %||% list())
    shown <- stats::setNames(vapply(outs, function(o) enc2utf8(o$displayed), ""),
                             vapply(outs, function(o) enc2utf8(o$label_as_displayed), ""))
    for (lab in names(cs$outputs)) expect_identical(shown[[enc2utf8(lab)]], cs$outputs[[lab]], info = paste(id, lab))
    if (!is.null(cs$message)) expect_true(enc2utf8(cs$message) %in% shown, info = id)
  }
})
