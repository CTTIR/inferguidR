# Curated G*Power fixture comparisons for group G1 (fast subset) ----------------------------
#
# A fast, offline subset of the emulator-capture comparisons of
# validation/compare/<method_id>/ (the full runs: validation/tools/
# compare_fixtures.R). Every case is a fixture record of
# validation/fixtures/emulator_capture/ or emulator_capture_supplement/g1/ (or
# registry-sync-and-comparator/), cited by its record id; the displayed
# strings below are copies of those records (G*Power 3.1.9.7 under Wine,
# evidence kind emulator_capture, decisions D22/D23), and the last test checks
# them against the fixture files in a development checkout. Every case listed
# here agrees in the full comparison; the records that disagree by design or
# through a pending core change are documented in the summaries and are not
# repeated here.
#
# Match rules (validation/tolerances.md, fixed before comparison): integers
# exact (MR-EXACT); displayed reals within the closed rounding interval of
# the last displayed digit plus 8 * eps * max(1, |x|) (MR-DISPLAY). G*Power
# refusals ('Input error') must be matched by a typed refusal without numbers
# (TD-0010).
#
# Display conventions (validation/compare/<method_id>/mapping.R): the t tests
# display the critical t on the side of delta, except in criterion mode,
# where G*Power shows the critical value of the mirrored test with |delta|;
# the generic z test shows the critical z on the side of mu in every mode;
# the effect-size drawers display |effect| except the matched-pairs mode
# 'From differences'.

ref_reset_registry()

g1_fixture_cases <- list(
  "GPW-T-MEAN1-APR-B001-O01--r01" = list(method = "t_means_one_sample", mode = "a_priori",
    inputs = c("Tail(s)" = "One", "Effect size d" = "0.625", "\u03b1 err prob" = "0.05", "Power (1-\u03b2 err prob)" = "0.95"),
    outputs = c("Noncentrality parameter \u03b4" = "3.4232660", "Critical t" = "1.6991270", "Df" = "29", "Total sample size" = "30", "Actual power" = "0.9551444")),
  "GPW-T-MEAN1-APR-B002-O01--r01" = list(method = "t_means_one_sample", mode = "a_priori",
    inputs = c("Tail(s)" = "Two", "Effect size d" = "0.1", "\u03b1 err prob" = "0.01", "Power (1-\u03b2 err prob)" = "0.90"),
    outputs = c("Noncentrality parameter \u03b4" = "3.8626416", "Critical t" = "2.5791308", "Df" = "1491", "Total sample size" = "1492", "Actual power" = "0.9001687")),
  "GPW-T-MEAN1-APR-B001-E01--r01" = list(method = "t_means_one_sample", mode = "a_priori",
    inputs = c("Tail(s)" = "One", "Effect size d" = "1.5", "\u03b1 err prob" = "0.05", "Power (1-\u03b2 err prob)" = "0.80"),
    outputs = c("Noncentrality parameter \u03b4" = "3.3541020", "Critical t" = "2.1318468", "Df" = "4", "Total sample size" = "5", "Actual power" = "0.8619466")),
  "SUP-G1-T-MEAN1-APR-A01--r01" = list(method = "t_means_one_sample", mode = "a_priori",
    inputs = c("Tail(s)" = "One", "Effect size d" = "10", "\u03b1 err prob" = "0.05", "Power (1-\u03b2 err prob)" = "0.80"),
    outputs = c("Noncentrality parameter \u03b4" = "14.1421356", "Critical t" = "6.3137515", "Df" = "1", "Total sample size" = "2", "Actual power" = "0.9730553")),
  "GPW-T-MEAN1-POS-B001-O01--r01" = list(method = "t_means_one_sample", mode = "post_hoc",
    inputs = c("Tail(s)" = "One", "Effect size d" = "0.5", "\u03b1 err prob" = "0.05", "Total sample size" = "34"),
    outputs = c("Noncentrality parameter \u03b4" = "2.9154759", "Critical t" = "1.6923603", "Df" = "33", "Power (1-\u03b2 err prob)" = "0.8868653")),
  "GPW-T-MEAN1-POS-B002-O01--r01" = list(method = "t_means_one_sample", mode = "post_hoc",
    inputs = c("Tail(s)" = "Two", "Effect size d" = "-0.25", "\u03b1 err prob" = "0.05", "Total sample size" = "25"),
    outputs = c("Noncentrality parameter \u03b4" = "-1.2500000", "Critical t" = "-2.0638986", "Df" = "24", "Power (1-\u03b2 err prob)" = "0.2245240")),
  "SUP-G1-T-MEAN1-POS-A01--r01" = list(method = "t_means_one_sample", mode = "post_hoc",
    inputs = c("Tail(s)" = "One", "Effect size d" = "0.5", "\u03b1 err prob" = "0.05", "Total sample size" = "3"),
    outputs = c("Noncentrality parameter \u03b4" = "0.8660254", "Critical t" = "2.9199856", "Df" = "2", "Power (1-\u03b2 err prob)" = "0.1512500")),
  "GPW-T-MEAN1-SEN-B001-O01--r01" = list(method = "t_means_one_sample", mode = "sensitivity",
    inputs = c("Tail(s)" = "One", "\u03b1 err prob" = "0.05", "Power (1-\u03b2 err prob)" = "0.80", "Total sample size" = "34"),
    outputs = c("Noncentrality parameter \u03b4" = "2.5390886", "Critical t" = "1.6923603", "Df" = "33", "Effect size d" = "0.4354501")),
  "GPW-T-MEAN1-SEN-B002-O01--r01" = list(method = "t_means_one_sample", mode = "sensitivity",
    inputs = c("Tail(s)" = "Two", "\u03b1 err prob" = "0.05", "Power (1-\u03b2 err prob)" = "0.80", "Total sample size" = "34"),
    outputs = c("Noncentrality parameter \u03b4" = "2.8864850", "Critical t" = "2.0345153", "Df" = "33", "Effect size d" = "0.4950281")),
  "GPW-T-MEAN1-CRI-B001-O01--r01" = list(method = "t_means_one_sample", mode = "criterion",
    inputs = c("Tail(s)" = "One", "Effect size d" = "0.5", "Power (1-\u03b2 err prob)" = "0.80", "Total sample size" = "34"),
    outputs = c("Noncentrality parameter \u03b4" = "2.9154759", "Critical t" = "2.0630050", "Df" = "33", "\u03b1 err prob" = "0.0235264")),
  "SUP-G1-T-MEAN1-CRI-A01--r01" = list(method = "t_means_one_sample", mode = "criterion",
    inputs = c("Tail(s)" = "One", "Effect size d" = "-0.5", "Power (1-\u03b2 err prob)" = "0.80", "Total sample size" = "34"),
    outputs = c("Noncentrality parameter \u03b4" = "-2.9154759", "Critical t" = "2.0630050", "Df" = "33", "\u03b1 err prob" = "0.0235264")),
  "SUP-G1-T-MEAN1-CRI-A02--r01" = list(method = "t_means_one_sample", mode = "criterion",
    inputs = c("Tail(s)" = "Two", "Effect size d" = "-0.5", "Power (1-\u03b2 err prob)" = "0.80", "Total sample size" = "34"),
    outputs = c("Noncentrality parameter \u03b4" = "-2.9154759", "Critical t" = "2.0630077", "Df" = "33", "\u03b1 err prob" = "0.0470526")),
  "GPW-T-MEAN1-COM-B001-O01--r01" = list(method = "t_means_one_sample", mode = "compromise",
    inputs = c("Tail(s)" = "One", "Effect size d" = "0.5", "\u03b2/\u03b1 ratio" = "1", "Total sample size" = "34"),
    outputs = c("Noncentrality parameter \u03b4" = "2.9154759", "Critical t" = "1.4686955", "Df" = "33", "\u03b1 err prob" = "0.0756921", "\u03b2 err prob" = "0.0756921", "Power (1-\u03b2 err prob)" = "0.9243079")),
  "GPW-T-MEAN1-COM-B002-O01--r01" = list(method = "t_means_one_sample", mode = "compromise",
    inputs = c("Tail(s)" = "Two", "Effect size d" = "0.5", "\u03b2/\u03b1 ratio" = "1", "Total sample size" = "34"),
    outputs = c("Noncentrality parameter \u03b4" = "2.9154759", "Critical t" = "1.6583345", "Df" = "33", "\u03b1 err prob" = "0.1067270", "\u03b2 err prob" = "0.1067270", "Power (1-\u03b2 err prob)" = "0.8932730")),
  "GPW-T-MEAN1-APR-B001-I03--r01" = list(method = "t_means_one_sample", mode = "a_priori",
    inputs = c("Tail(s)" = "One", "Effect size d" = "0.5", "\u03b1 err prob" = "0.05", "Power (1-\u03b2 err prob)" = "0.04"),
    outputs = c("Noncentrality parameter \u03b4" = "?", "Critical t" = "?", "Df" = "?", "Total sample size" = "?", "Actual power" = "?"),
    message = "Power must be greater than or equal to alpha."),
  "GPW-T-MEAN1-SEN-B001-I02--r01" = list(method = "t_means_one_sample", mode = "sensitivity",
    inputs = c("Tail(s)" = "One", "\u03b1 err prob" = "0.05", "Power (1-\u03b2 err prob)" = "0.04", "Total sample size" = "34"),
    outputs = c("Noncentrality parameter \u03b4" = "?", "Critical t" = "?", "Df" = "?", "Effect size d" = "?"),
    message = "Power must be greater than or equal to alpha."),
  "GPW-T-MEAN1-COM-B001-I01--r01" = list(method = "t_means_one_sample", mode = "compromise",
    inputs = c("Tail(s)" = "One", "Effect size d" = "0.5", "\u03b2/\u03b1 ratio" = "0", "Total sample size" = "34"),
    outputs = c("Noncentrality parameter \u03b4" = "?", "Critical t" = "?", "Df" = "?", "\u03b1 err prob" = "?", "\u03b2 err prob" = "?", "Power (1-\u03b2 err prob)" = "?"),
    message = "\u03b2/\u03b1 ratio is outside the admissible range (0.01, 100)"),
  "GPW-UT-DRAW-UTL-U044-UO01--r01" = list(method = "t_means_one_sample", mode = "effect_size_calculator",
    inputs = c("Tail(s)" = "One", "Effect size d" = "0.5", "\u03b1 err prob" = "0.05", "Power (1-\u03b2 err prob)" = "0.95", "Mean H0" = "10", "Mean H1" = "15", "SD \u03c3" = "8", "drawer mode" = "no selector matched (drawer mode '(mode not printed; select the mode named by the branch)' not found among the live combo items, radios and group boxes); fields used as displayed"),
    outputs = c("Effect size d" = "0.625")),
  "SUP-G1-DRAW-MEAN1-UD01--r01" = list(method = "t_means_one_sample", mode = "effect_size_calculator",
    inputs = c("Tail(s)" = "One", "Effect size d" = "0.5", "\u03b1 err prob" = "0.05", "Power (1-\u03b2 err prob)" = "0.95", "Mean H0" = "15", "Mean H1" = "10", "SD \u03c3" = "8", "drawer mode" = "no selector matched (drawer mode '(mode not printed; select the mode named by the branch)' not found among the live combo items, radios and group boxes); fields used as displayed"),
    outputs = c("Effect size d" = "0.625")),
  "GPW-T-PAIRED-APR-B001-O01--r01" = list(method = "t_means_matched_pairs", mode = "a_priori",
    inputs = c("Tail(s)" = "One", "Effect size dz" = "0.5", "\u03b1 err prob" = "0.05", "Power (1-\u03b2 err prob)" = "0.80"),
    outputs = c("Noncentrality parameter \u03b4" = "2.5980762", "Critical t" = "1.7056179", "Df" = "26", "Total sample size" = "27", "Actual power" = "0.8118316")),
  "GPW-T-PAIRED-APR-B002-D01--r01" = list(method = "t_means_matched_pairs", mode = "a_priori",
    inputs = c("Tail(s)" = "Two", "Effect size dz" = "0.02", "\u03b1 err prob" = "0.001", "Power (1-\u03b2 err prob)" = "0.99"),
    outputs = c("Noncentrality parameter \u03b4" = "5.6170811", "Critical t" = "3.2906501", "Df" = "78878", "Total sample size" = "78879", "Actual power" = "0.9900004")),
  "SUP-G1-T-PAIRED-APR-A01--r01" = list(method = "t_means_matched_pairs", mode = "a_priori",
    inputs = c("Tail(s)" = "One", "Effect size dz" = "10", "\u03b1 err prob" = "0.05", "Power (1-\u03b2 err prob)" = "0.80"),
    outputs = c("Noncentrality parameter \u03b4" = "14.1421356", "Critical t" = "6.3137515", "Df" = "1", "Total sample size" = "2", "Actual power" = "0.9730553")),
  "GPW-T-PAIRED-POS-B002-O01--r01" = list(method = "t_means_matched_pairs", mode = "post_hoc",
    inputs = c("Tail(s)" = "Two", "Effect size dz" = "0.421637", "\u03b1 err prob" = "0.05", "Total sample size" = "50"),
    outputs = c("Noncentrality parameter \u03b4" = "2.9814238", "Critical t" = "2.0095752", "Df" = "49", "Power (1-\u03b2 err prob)" = "0.8321145")),
  "GPW-T-PAIRED-SEN-B002-D01--r01" = list(method = "t_means_matched_pairs", mode = "sensitivity",
    inputs = c("Tail(s)" = "Two", "\u03b1 err prob" = "0.001", "Power (1-\u03b2 err prob)" = "0.99", "Total sample size" = "100000"),
    outputs = c("Noncentrality parameter \u03b4" = "5.6170267", "Critical t" = "3.2906240", "Df" = "99999", "Effect size dz" = "0.0177626")),
  "GPW-T-PAIRED-CRI-B002-O01--r01" = list(method = "t_means_matched_pairs", mode = "criterion",
    inputs = c("Tail(s)" = "Two", "Effect size dz" = "0.5", "Power (1-\u03b2 err prob)" = "0.80", "Total sample size" = "34"),
    outputs = c("Noncentrality parameter \u03b4" = "2.9154759", "Critical t" = "2.0630077", "Df" = "33", "\u03b1 err prob" = "0.0470526")),
  "SUP-G1-T-PAIRED-CRI-A01--r01" = list(method = "t_means_matched_pairs", mode = "criterion",
    inputs = c("Tail(s)" = "One", "Effect size dz" = "-0.5", "Power (1-\u03b2 err prob)" = "0.80", "Total sample size" = "34"),
    outputs = c("Noncentrality parameter \u03b4" = "-2.9154759", "Critical t" = "2.0630050", "Df" = "33", "\u03b1 err prob" = "0.0235264")),
  "GPW-T-PAIRED-COM-B001-O01--r01" = list(method = "t_means_matched_pairs", mode = "compromise",
    inputs = c("Tail(s)" = "One", "Effect size dz" = "0.5", "\u03b2/\u03b1 ratio" = "1", "Total sample size" = "34"),
    outputs = c("Noncentrality parameter \u03b4" = "2.9154759", "Critical t" = "1.4686955", "Df" = "33", "\u03b1 err prob" = "0.0756921", "\u03b2 err prob" = "0.0756921", "Power (1-\u03b2 err prob)" = "0.9243079")),
  "GPW-UT-DRAW-UTL-U045-UE01--r01" = list(method = "t_means_matched_pairs", mode = "effect_size_calculator",
    inputs = c("Tail(s)" = "One", "Effect size dz" = "0.5", "\u03b1 err prob" = "0.05", "Power (1-\u03b2 err prob)" = "0.95", "Mean of difference" = "-5", "SD of difference" = "0.1", "drawer mode" = "group 'From differences'"),
    outputs = c("Effect size dz" = "-50")),
  "GPW-UT-DRAW-UTL-U046-UO01--r01" = list(method = "t_means_matched_pairs", mode = "effect_size_calculator",
    inputs = c("Tail(s)" = "One", "Effect size dz" = "0.5", "\u03b1 err prob" = "0.05", "Power (1-\u03b2 err prob)" = "0.95", "Mean group 1" = "0", "Mean group 2" = "0.4", "SD group 1" = "1", "SD group 2" = "1", "Correlation between groups" = "0.55", "drawer mode" = "group 'From group parameters'"),
    outputs = c("Effect size dz" = "0.421637")),
  "SUP-G1-DRAW-PAIRED-UD01--r01" = list(method = "t_means_matched_pairs", mode = "effect_size_calculator",
    inputs = c("Tail(s)" = "One", "Effect size dz" = "0.5", "\u03b1 err prob" = "0.05", "Power (1-\u03b2 err prob)" = "0.95", "Mean group 1" = "0.4", "Mean group 2" = "0", "SD group 1" = "1", "SD group 2" = "1", "Correlation between groups" = "0.55", "drawer mode" = "group 'From group parameters'"),
    outputs = c("Effect size dz" = "0.421637")),
  "SUP-RSC-T-MEAN2-APR-A01--r01" = list(method = "t_means_two_independent", mode = "a_priori",
    inputs = c("Tail(s)" = "One", "Effect size d" = "0.5", "\u03b1 err prob" = "0.05", "Power (1-\u03b2 err prob)" = "0.80", "Allocation ratio N2/N1" = "3"),
    outputs = c("Noncentrality parameter \u03b4" = "2.5185876", "Critical t" = "1.6564793", "Df" = "132", "Sample size group 1" = "34", "Sample size group 2" = "100", "Total sample size" = "134", "Actual power" = "0.8053237")),
  "SUP-RSC-T-MEAN2-APR-A02--r01" = list(method = "t_means_two_independent", mode = "a_priori",
    inputs = c("Tail(s)" = "Two", "Effect size d" = "0.5", "\u03b1 err prob" = "0.05", "Power (1-\u03b2 err prob)" = "0.80", "Allocation ratio N2/N1" = "3"),
    outputs = c("Noncentrality parameter \u03b4" = "2.8338811", "Critical t" = "1.9741852", "Df" = "168", "Sample size group 1" = "43", "Sample size group 2" = "127", "Total sample size" = "170", "Actual power" = "0.8044668")),
  "SUP-RSC-T-MEAN2-APR-A03--r01" = list(method = "t_means_two_independent", mode = "a_priori",
    inputs = c("Tail(s)" = "One", "Effect size d" = "0.8", "\u03b1 err prob" = "0.05", "Power (1-\u03b2 err prob)" = "0.80", "Allocation ratio N2/N1" = "0.2"),
    outputs = c("Noncentrality parameter \u03b4" = "2.5298221", "Critical t" = "1.6669145", "Df" = "70", "Sample size group 1" = "60", "Sample size group 2" = "12", "Total sample size" = "72", "Actual power" = "0.8052141")),
  "SUP-RSC-T-MEAN2-APR-A05--r01" = list(method = "t_means_two_independent", mode = "a_priori",
    inputs = c("Tail(s)" = "Two", "Effect size d" = "1.5", "\u03b1 err prob" = "0.05", "Power (1-\u03b2 err prob)" = "0.95", "Allocation ratio N2/N1" = "2.5"),
    outputs = c("Noncentrality parameter \u03b4" = "3.8150606", "Critical t" = "2.0422725", "Df" = "30", "Sample size group 1" = "9", "Sample size group 2" = "23", "Total sample size" = "32", "Actual power" = "0.9582163")),
  "GPW-T-MEAN2-APR-B001-D01--r01" = list(method = "t_means_two_independent", mode = "a_priori",
    inputs = c("Tail(s)" = "One", "Effect size d" = "0.02", "\u03b1 err prob" = "0.001", "Power (1-\u03b2 err prob)" = "0.99", "Allocation ratio N2/N1" = "0.1"),
    outputs = c("Noncentrality parameter \u03b4" = "5.4166062", "Critical t" = "3.0902415", "Df" = "887518", "Sample size group 1" = "806836", "Sample size group 2" = "80684", "Total sample size" = "887520", "Actual power" = "0.9900003")),
  "GPW-T-MEAN2-APR-B002-E01--r01" = list(method = "t_means_two_independent", mode = "a_priori",
    inputs = c("Tail(s)" = "Two", "Effect size d" = "1.5", "\u03b1 err prob" = "0.05", "Power (1-\u03b2 err prob)" = "0.80", "Allocation ratio N2/N1" = "10"),
    outputs = c("Noncentrality parameter \u03b4" = "2.8665993", "Critical t" = "2.0153676", "Df" = "44", "Sample size group 1" = "4", "Sample size group 2" = "42", "Total sample size" = "46", "Actual power" = "0.8005431")),
  "GPW-T-MEAN2-POS-B001-E01--r01" = list(method = "t_means_two_independent", mode = "post_hoc",
    inputs = c("Tail(s)" = "One", "Effect size d" = "-0.5", "\u03b1 err prob" = "0.05", "Sample size group 1" = "2", "Sample size group 2" = "20"),
    outputs = c("Noncentrality parameter \u03b4" = "-0.6741999", "Critical t" = "-1.7247182", "Df" = "20", "Power (1-\u03b2 err prob)" = "0.1602865")),
  "GPW-T-MEAN2-SEN-B002-E01--r01" = list(method = "t_means_two_independent", mode = "sensitivity",
    inputs = c("Tail(s)" = "Two", "\u03b1 err prob" = "0.05", "Power (1-\u03b2 err prob)" = "0.95", "Sample size group 1" = "2", "Sample size group 2" = "20"),
    outputs = c("Noncentrality parameter \u03b4" = "3.7924567", "Critical t" = "2.0859634", "Df" = "20", "Effect size d" = "2.8125611")),
  "GPW-T-MEAN2-CRI-B001-E01--r01" = list(method = "t_means_two_independent", mode = "criterion",
    inputs = c("Tail(s)" = "One", "Effect size d" = "1.5", "Power (1-\u03b2 err prob)" = "0.95", "Sample size group 1" = "2", "Sample size group 2" = "20"),
    outputs = c("Noncentrality parameter \u03b4" = "2.0225996", "Critical t" = "0.3795277", "Df" = "20", "\u03b1 err prob" = "0.3541480")),
  "GPW-T-MEAN2-COM-B002-E01--r01" = list(method = "t_means_two_independent", mode = "compromise",
    inputs = c("Tail(s)" = "Two", "Effect size d" = "1.5", "\u03b2/\u03b1 ratio" = "4", "Sample size group 1" = "2", "Sample size group 2" = "20"),
    outputs = c("Noncentrality parameter \u03b4" = "2.0225996", "Critical t" = "1.7442762", "Df" = "20", "\u03b1 err prob" = "0.0964613", "\u03b2 err prob" = "0.3858453", "Power (1-\u03b2 err prob)" = "0.6141547")),
  "GPW-T-MEAN2-APR-B001-I03--r01" = list(method = "t_means_two_independent", mode = "a_priori",
    inputs = c("Tail(s)" = "One", "Effect size d" = "0.5", "\u03b1 err prob" = "0.05", "Power (1-\u03b2 err prob)" = "0.04", "Allocation ratio N2/N1" = "1"),
    outputs = c("Noncentrality parameter \u03b4" = "?", "Critical t" = "?", "Df" = "?", "Sample size group 1" = "?", "Sample size group 2" = "?", "Total sample size" = "?", "Actual power" = "?"),
    message = "Power must be greater than or equal to alpha."),
  "SUP-G1-DRAW-MEAN2-UD01--r01" = list(method = "t_means_two_independent", mode = "effect_size_calculator",
    inputs = c("Tail(s)" = "One", "Effect size d" = "0.5", "\u03b1 err prob" = "0.05", "Power (1-\u03b2 err prob)" = "0.95", "Allocation ratio N2/N1" = "1", "Mean group 1" = "1", "Mean group 2" = "0", "SD \u03c3 within each group" = "2", "drawer mode" = "group 'n1 != n2'"),
    outputs = c("Effect size d" = "0.5")),
  "SUP-G1-DRAW-MEAN2-UD02--r01" = list(method = "t_means_two_independent", mode = "effect_size_calculator",
    inputs = c("Tail(s)" = "One", "Effect size d" = "0.5", "\u03b1 err prob" = "0.05", "Power (1-\u03b2 err prob)" = "0.95", "Allocation ratio N2/N1" = "1", "Mean group 1" = "1", "Mean group 2" = "0", "SD \u03c3 group 1" = "0.5", "SD \u03c3 group 2" = "0.5", "drawer mode" = "group 'n1 = n2'"),
    outputs = c("Effect size d" = "2")),
  "GPW-T-GENT-POS-B001-O01--r01" = list(method = "t_generic", mode = "post_hoc",
    inputs = c("Tail(s)" = "One", "Noncentrality parameter \u03b4" = "2.5", "\u03b1 err prob" = "0.05", "Df" = "30"),
    outputs = c("Critical t" = "1.6972609", "\u03b2 err prob" = "0.2123992", "Power (1-\u03b2 err prob)" = "0.7876008")),
  "GPW-T-GENT-POS-B002-E01--r01" = list(method = "t_generic", mode = "post_hoc",
    inputs = c("Tail(s)" = "Two", "Noncentrality parameter \u03b4" = "-1", "\u03b1 err prob" = "0.05", "Df" = "1"),
    outputs = c("Critical t" = "-12.7062047", "\u03b2 err prob" = "0.9270111", "Power (1-\u03b2 err prob)" = "0.0729889")),
  "GPW-T-GENT-POS-B001-D01--r01" = list(method = "t_generic", mode = "post_hoc",
    inputs = c("Tail(s)" = "One", "Noncentrality parameter \u03b4" = "40", "\u03b1 err prob" = "0.0001", "Df" = "100000"),
    outputs = c("Critical t" = "3.7191544", "\u03b2 err prob" = "0", "Power (1-\u03b2 err prob)" = "1.0000000")),
  "SUP-G1-T-GENT-POS-A01--r01" = list(method = "t_generic", mode = "post_hoc",
    inputs = c("Tail(s)" = "One", "Noncentrality parameter \u03b4" = "2.5", "\u03b1 err prob" = "0.05", "Df" = "0.5"),
    outputs = c("Critical t" = "41.1360001", "\u03b2 err prob" = "0.8125006", "Power (1-\u03b2 err prob)" = "0.1874994")),
  "GPW-T-GENT-SEN-B001-E01--r01" = list(method = "t_generic", mode = "sensitivity",
    inputs = c("Tail(s)" = "One", "\u03b1 err prob" = "0.05", "Power (1-\u03b2 err prob)" = "0.95", "Df" = "1"),
    outputs = c("Critical t" = "6.3137515", "Noncentrality parameter \u03b4" = "12.5289781")),
  "GPW-T-GENT-SEN-B002-O01--r01" = list(method = "t_generic", mode = "sensitivity",
    inputs = c("Tail(s)" = "Two", "\u03b1 err prob" = "0.05", "Power (1-\u03b2 err prob)" = "0.80", "Df" = "30"),
    outputs = c("Critical t" = "2.0422725", "Noncentrality parameter \u03b4" = "2.8953583")),
  "GPW-T-GENT-CRI-B001-E01--r01" = list(method = "t_generic", mode = "criterion",
    inputs = c("Tail(s)" = "One", "Noncentrality parameter \u03b4" = "-1", "Power (1-\u03b2 err prob)" = "0.95", "Df" = "1"),
    outputs = c("Critical t" = "-1.1931118", "\u03b1 err prob" = "0.7779563")),
  "GPW-T-GENT-CRI-B002-E01--r01" = list(method = "t_generic", mode = "criterion",
    inputs = c("Tail(s)" = "Two", "Noncentrality parameter \u03b4" = "-1", "Power (1-\u03b2 err prob)" = "0.95", "Df" = "1"),
    outputs = c("Critical t" = "0.1294950", "\u03b1 err prob" = "0.9180171")),
  "SUP-G1-T-GENT-CRI-A01--r01" = list(method = "t_generic", mode = "criterion",
    inputs = c("Tail(s)" = "One", "Noncentrality parameter \u03b4" = "-2.5", "Power (1-\u03b2 err prob)" = "0.80", "Df" = "30"),
    outputs = c("Critical t" = "1.6533561", "\u03b1 err prob" = "0.0543434")),
  "GPW-T-GENT-COM-B001-E01--r01" = list(method = "t_generic", mode = "compromise",
    inputs = c("Tail(s)" = "One", "Noncentrality parameter \u03b4" = "-1", "\u03b2/\u03b1 ratio" = "4", "Df" = "1"),
    outputs = c("Critical t" = "-1.9207405", "\u03b1 err prob" = "0.1527942", "\u03b2 err prob" = "0.6111767", "Power (1-\u03b2 err prob)" = "0.3888233")),
  "GPW-T-GENT-COM-B002-O01--r01" = list(method = "t_generic", mode = "compromise",
    inputs = c("Tail(s)" = "Two", "Noncentrality parameter \u03b4" = "2.5", "\u03b2/\u03b1 ratio" = "1", "Df" = "30"),
    outputs = c("Critical t" = "1.4682587", "\u03b1 err prob" = "0.1524418", "\u03b2 err prob" = "0.1524418", "Power (1-\u03b2 err prob)" = "0.8475582")),
  "GPW-T-GENT-POS-B001-I02--r01" = list(method = "t_generic", mode = "post_hoc",
    inputs = c("Tail(s)" = "One", "Noncentrality parameter \u03b4" = "2.5", "\u03b1 err prob" = "0.05", "Df" = "0"),
    outputs = c("Critical t" = "?", "\u03b2 err prob" = "?", "Power (1-\u03b2 err prob)" = "?"),
    message = "Degrees of Freedom must be greater 0.01."),
  "GPW-Z-GENZ-POS-B001-E01--r01" = list(method = "z_generic", mode = "post_hoc",
    inputs = c("Tail(s)" = "One", "Noncentrality parameter \u03bc" = "-0.5", "Noncentral dist. SD \u03c3" = "0.2", "\u03b1 err prob" = "0.05"),
    outputs = c("Critical z" = "-1.6448536", "\u03b2 err prob" = "1.0000000", "Power (1-\u03b2 err prob)" = "5.194029e-09")),
  "GPW-Z-GENZ-POS-B002-E01--r01" = list(method = "z_generic", mode = "post_hoc",
    inputs = c("Tail(s)" = "Two", "Noncentrality parameter \u03bc" = "-0.5", "Noncentral dist. SD \u03c3" = "0.2", "\u03b1 err prob" = "0.05"),
    outputs = c("Critical z" = "-1.9599640", "\u03b2 err prob" = "1.0000000", "Power (1-\u03b2 err prob)" = "1.441069e-13")),
  "GPW-Z-GENZ-POS-B002-D01--r01" = list(method = "z_generic", mode = "post_hoc",
    inputs = c("Tail(s)" = "Two", "Noncentrality parameter \u03bc" = "0.001", "Noncentral dist. SD \u03c3" = "5", "\u03b1 err prob" = "0.0001"),
    outputs = c("Critical z" = "3.8905919", "\u03b2 err prob" = "0.5635008", "Power (1-\u03b2 err prob)" = "0.4364992")),
  "GPW-Z-GENZ-SEN-B001-D01--r01" = list(method = "z_generic", mode = "sensitivity",
    inputs = c("Tail(s)" = "One", "Noncentral dist. SD \u03c3" = "5", "\u03b1 err prob" = "0.001", "Power (1-\u03b2 err prob)" = "0.99"),
    outputs = c("Critical z" = "3.0902323", "Noncentrality parameter \u03bc" = "14.7219717")),
  "GPW-Z-GENZ-SEN-B002-E01--r01" = list(method = "z_generic", mode = "sensitivity",
    inputs = c("Tail(s)" = "Two", "Noncentral dist. SD \u03c3" = "0.2", "\u03b1 err prob" = "0.05", "Power (1-\u03b2 err prob)" = "0.95"),
    outputs = c("Critical z" = "1.9599640", "Noncentrality parameter \u03bc" = "2.2889347")),
  "GPW-Z-GENZ-CRI-B001-E01--r01" = list(method = "z_generic", mode = "criterion",
    inputs = c("Tail(s)" = "One", "Noncentrality parameter \u03bc" = "6", "Noncentral dist. SD \u03c3" = "0.2", "Power (1-\u03b2 err prob)" = "0.95"),
    outputs = c("Critical z" = "5.6710293", "\u03b1 err prob" = "7.097107e-09")),
  "GPW-Z-GENZ-CRI-B002-D01--r01" = list(method = "z_generic", mode = "criterion",
    inputs = c("Tail(s)" = "Two", "Noncentrality parameter \u03bc" = "0.001", "Noncentral dist. SD \u03c3" = "5", "Power (1-\u03b2 err prob)" = "0.9999"),
    outputs = c("Critical z" = "0.0006266571", "\u03b1 err prob" = "0.9995000")),
  "SUP-G1-Z-GENZ-CRI-A01--r01" = list(method = "z_generic", mode = "criterion",
    inputs = c("Tail(s)" = "One", "Noncentrality parameter \u03bc" = "-2.5", "Noncentral dist. SD \u03c3" = "1", "Power (1-\u03b2 err prob)" = "0.80"),
    outputs = c("Critical z" = "-1.6583788", "\u03b1 err prob" = "0.0486205")),
  "SUP-G1-Z-GENZ-CRI-A02--r01" = list(method = "z_generic", mode = "criterion",
    inputs = c("Tail(s)" = "Two", "Noncentrality parameter \u03bc" = "-2.5", "Noncentral dist. SD \u03c3" = "1", "Power (1-\u03b2 err prob)" = "0.80"),
    outputs = c("Critical z" = "-1.6584360", "\u03b1 err prob" = "0.0972295")),
  "GPW-Z-GENZ-COM-B001-E01--r01" = list(method = "z_generic", mode = "compromise",
    inputs = c("Tail(s)" = "One", "Noncentrality parameter \u03bc" = "6", "Noncentral dist. SD \u03c3" = "0.2", "\u03b2/\u03b1 ratio" = "4"),
    outputs = c("Critical z" = "5.0453133", "\u03b1 err prob" = "2.263894e-07", "\u03b2 err prob" = "9.055575e-07", "Power (1-\u03b2 err prob)" = "0.9999991")),
  "GPW-Z-GENZ-COM-B002-D01--r01" = list(method = "z_generic", mode = "compromise",
    inputs = c("Tail(s)" = "Two", "Noncentrality parameter \u03bc" = "0.001", "Noncentral dist. SD \u03c3" = "5", "\u03b2/\u03b1 ratio" = "0.01"),
    outputs = c("Critical z" = "0.0596847", "\u03b1 err prob" = "0.9524068", "\u03b2 err prob" = "0.009524068", "Power (1-\u03b2 err prob)" = "0.9904759")),
  "GPW-Z-GENZ-POS-B001-I02--r01" = list(method = "z_generic", mode = "post_hoc",
    inputs = c("Tail(s)" = "One", "Noncentrality parameter \u03bc" = "2.5", "Noncentral dist. SD \u03c3" = "0", "\u03b1 err prob" = "0.05"),
    outputs = c("Critical z" = "?", "\u03b2 err prob" = "?", "Power (1-\u03b2 err prob)" = "?"),
    message = "Standard error \u03c3 must be greater or equal 1e-08."),
  "GPW-F-GENF-POS-B001-O01--r01" = list(method = "f_generic", mode = "post_hoc",
    inputs = c("Noncentrality parameter \u03bb" = "12", "\u03b1 err prob" = "0.05", "Numerator df" = "3", "Denominator df" = "60"),
    outputs = c("Critical F" = "2.7580783", "\u03b2 err prob" = "0.1867858", "Power (1-\u03b2 err prob)" = "0.8132142")),
  "GPW-F-GENF-POS-B001-D01--r01" = list(method = "f_generic", mode = "post_hoc",
    inputs = c("Noncentrality parameter \u03bb" = "400", "\u03b1 err prob" = "0.0001", "Numerator df" = "50", "Denominator df" = "100000"),
    outputs = c("Critical F" = "1.9198353", "\u03b2 err prob" = "0", "Power (1-\u03b2 err prob)" = "1.0000000")),
  "GPW-F-GENF-POS-B001-E01--r01" = list(method = "f_generic", mode = "post_hoc",
    inputs = c("Noncentrality parameter \u03bb" = "0.5", "\u03b1 err prob" = "0.05", "Numerator df" = "1", "Denominator df" = "1"),
    outputs = c("Critical F" = "161.4476", "\u03b2 err prob" = "0.9380514", "Power (1-\u03b2 err prob)" = "0.0619486")),
  "SUP-G1-F-GENF-POS-A01--r01" = list(method = "f_generic", mode = "post_hoc",
    inputs = c("Noncentrality parameter \u03bb" = "12", "\u03b1 err prob" = "0.05", "Numerator df" = "2.5", "Denominator df" = "7.5"),
    outputs = c("Critical F" = "4.3600974", "\u03b2 err prob" = "0.3479669", "Power (1-\u03b2 err prob)" = "0.6520331")),
  "GPW-F-GENF-SEN-B001-O01--r01" = list(method = "f_generic", mode = "sensitivity",
    inputs = c("\u03b1 err prob" = "0.05", "Power (1-\u03b2 err prob)" = "0.80", "Numerator df" = "3", "Denominator df" = "60"),
    outputs = c("Critical F" = "2.7580783", "Noncentrality parameter \u03bb" = "11.6382121")),
  "GPW-F-GENF-SEN-B001-D01--r01" = list(method = "f_generic", mode = "sensitivity",
    inputs = c("\u03b1 err prob" = "0.001", "Power (1-\u03b2 err prob)" = "0.99", "Numerator df" = "50", "Denominator df" = "100000"),
    outputs = c("Critical F" = "1.7335513", "Noncentrality parameter \u03bb" = "80.2863927")),
  "SUP-G1-F-GENF-SEN-A01--r01" = list(method = "f_generic", mode = "sensitivity",
    inputs = c("\u03b1 err prob" = "0.05", "Power (1-\u03b2 err prob)" = "0.56", "Numerator df" = "1", "Denominator df" = "1"),
    outputs = c("Critical F" = "161.4476", "Noncentrality parameter \u03bb" = "96.8646614")),
  "GPW-F-GENF-CRI-B001-O01--r01" = list(method = "f_generic", mode = "criterion",
    inputs = c("Noncentrality parameter \u03bb" = "12", "Power (1-\u03b2 err prob)" = "0.80", "Numerator df" = "3", "Denominator df" = "60"),
    outputs = c("Critical F" = "2.8482527", "\u03b1 err prob" = "0.0449105")),
  "GPW-F-GENF-CRI-B001-E01--r01" = list(method = "f_generic", mode = "criterion",
    inputs = c("Noncentrality parameter \u03bb" = "0.5", "Power (1-\u03b2 err prob)" = "0.95", "Numerator df" = "1", "Denominator df" = "1"),
    outputs = c("Critical F" = "0.0102048", "\u03b1 err prob" = "0.9359069")),
  "GPW-F-GENF-COM-B001-O01--r01" = list(method = "f_generic", mode = "compromise",
    inputs = c("Noncentrality parameter \u03bb" = "12", "\u03b2/\u03b1 ratio" = "1", "Numerator df" = "3", "Denominator df" = "60"),
    outputs = c("Critical F" = "2.1342916", "\u03b1 err prob" = "0.1052850", "\u03b2 err prob" = "0.1052850", "Power (1-\u03b2 err prob)" = "0.8947150")),
  "GPW-F-GENF-COM-B001-E01--r01" = list(method = "f_generic", mode = "compromise",
    inputs = c("Noncentrality parameter \u03bb" = "0.5", "\u03b2/\u03b1 ratio" = "4", "Numerator df" = "1", "Denominator df" = "1"),
    outputs = c("Critical F" = "10.4006582", "\u03b1 err prob" = "0.1914162", "\u03b2 err prob" = "0.7656650", "Power (1-\u03b2 err prob)" = "0.2343350")),
  "GPW-F-GENF-POS-B001-I02--r01" = list(method = "f_generic", mode = "post_hoc",
    inputs = c("Noncentrality parameter \u03bb" = "-1", "\u03b1 err prob" = "0.05", "Numerator df" = "3", "Denominator df" = "60"),
    outputs = c("Critical F" = "?", "\u03b2 err prob" = "?", "Power (1-\u03b2 err prob)" = "?"),
    message = "Noncentrality parameter \u03bb must be non-negative."),
  "GPW-CHI2-GENCHI-POS-B001-O01--r01" = list(method = "chi2_generic", mode = "post_hoc",
    inputs = c("Noncentrality parameter \u03bb" = "10", "\u03b1 err prob" = "0.05", "Df" = "5"),
    outputs = c("Critical \u03c7\u00b2" = "11.0704977", "\u03b2 err prob" = "0.3225611", "Power (1-\u03b2 err prob)" = "0.6774389")),
  "GPW-CHI2-GENCHI-POS-B001-D01--r01" = list(method = "chi2_generic", mode = "post_hoc",
    inputs = c("Noncentrality parameter \u03bb" = "400", "\u03b1 err prob" = "0.0001", "Df" = "1000"),
    outputs = c("Critical \u03c7\u00b2" = "1174.933", "\u03b2 err prob" = "3.935676e-05", "Power (1-\u03b2 err prob)" = "0.9999606")),
  "SUP-G1-CHI2-GENCHI-POS-A01--r01" = list(method = "chi2_generic", mode = "post_hoc",
    inputs = c("Noncentrality parameter \u03bb" = "10", "\u03b1 err prob" = "0.05", "Df" = "0.5"),
    outputs = c("Critical \u03c7\u00b2" = "2.4202323", "\u03b2 err prob" = "0.0678334", "Power (1-\u03b2 err prob)" = "0.9321666")),
  "GPW-CHI2-GENCHI-SEN-B001-E01--r01" = list(method = "chi2_generic", mode = "sensitivity",
    inputs = c("\u03b1 err prob" = "0.05", "Power (1-\u03b2 err prob)" = "0.95", "Df" = "1"),
    outputs = c("Critical \u03c7\u00b2" = "3.8414588", "Noncentrality parameter \u03bb" = "12.9947091")),
  "SUP-G1-CHI2-GENCHI-SEN-A02--r01" = list(method = "chi2_generic", mode = "sensitivity",
    inputs = c("\u03b1 err prob" = "0.05", "Power (1-\u03b2 err prob)" = "0.69", "Df" = "1000"),
    outputs = c("Critical \u03c7\u00b2" = "1074.679", "Noncentrality parameter \u03bb" = "99.4965539")),
  "GPW-CHI2-GENCHI-CRI-B001-D01--r01" = list(method = "chi2_generic", mode = "criterion",
    inputs = c("Noncentrality parameter \u03bb" = "400", "Power (1-\u03b2 err prob)" = "0.9999", "Df" = "1000"),
    outputs = c("Critical \u03c7\u00b2" = "1187.274", "\u03b1 err prob" = "3.652135e-05")),
  "GPW-CHI2-GENCHI-CRI-B001-E01--r01" = list(method = "chi2_generic", mode = "criterion",
    inputs = c("Noncentrality parameter \u03bb" = "0.5", "Power (1-\u03b2 err prob)" = "0.95", "Df" = "1"),
    outputs = c("Critical \u03c7\u00b2" = "0.006481512", "\u03b1 err prob" = "0.9358334")),
  "GPW-CHI2-GENCHI-COM-B001-D01--r01" = list(method = "chi2_generic", mode = "compromise",
    inputs = c("Noncentrality parameter \u03bb" = "400", "\u03b2/\u03b1 ratio" = "0.01", "Df" = "1000"),
    outputs = c("Critical \u03c7\u00b2" = "1151.666", "\u03b1 err prob" = "0.0005803717", "\u03b2 err prob" = "5.803717e-06", "Power (1-\u03b2 err prob)" = "0.9999942")),
  "GPW-CHI2-GENCHI-COM-B001-O01--r01" = list(method = "chi2_generic", mode = "compromise",
    inputs = c("Noncentrality parameter \u03bb" = "10", "\u03b2/\u03b1 ratio" = "1", "Df" = "5"),
    outputs = c("Critical \u03c7\u00b2" = "8.0107153", "\u03b1 err prob" = "0.1556461", "\u03b2 err prob" = "0.1556461", "Power (1-\u03b2 err prob)" = "0.8443539")),
  "GPW-CHI2-GENCHI-SEN-B001-I02--r01" = list(method = "chi2_generic", mode = "sensitivity",
    inputs = c("\u03b1 err prob" = "0.05", "Power (1-\u03b2 err prob)" = "0.04", "Df" = "5"),
    outputs = c("Critical \u03c7\u00b2" = "?", "Noncentrality parameter \u03bb" = "?"),
    message = "Power must be greater than or equal to alpha.")
)

g1_input_map <- list(
  t_means_one_sample = c("Effect size d" = "d", "\u03b1 err prob" = "alpha", "Power (1-\u03b2 err prob)" = "power",
                         "\u03b2/\u03b1 ratio" = "error_ratio", "Total sample size" = "n_total"),
  t_means_matched_pairs = c("Effect size dz" = "dz", "\u03b1 err prob" = "alpha", "Power (1-\u03b2 err prob)" = "power",
                            "\u03b2/\u03b1 ratio" = "error_ratio", "Total sample size" = "n_total"),
  t_means_two_independent = c("Effect size d" = "d", "\u03b1 err prob" = "alpha", "Power (1-\u03b2 err prob)" = "power",
                              "\u03b2/\u03b1 ratio" = "error_ratio", "Sample size group 1" = "n1",
                              "Sample size group 2" = "n2", "Allocation ratio N2/N1" = "allocation_ratio"),
  t_generic = c("Noncentrality parameter \u03b4" = "delta", "\u03b1 err prob" = "alpha",
                "Power (1-\u03b2 err prob)" = "power", "\u03b2/\u03b1 ratio" = "error_ratio", "Df" = "df"),
  z_generic = c("Noncentrality parameter \u03bc" = "mu", "Noncentral dist. SD \u03c3" = "sigma",
                "\u03b1 err prob" = "alpha", "Power (1-\u03b2 err prob)" = "power", "\u03b2/\u03b1 ratio" = "error_ratio"),
  f_generic = c("Noncentrality parameter \u03bb" = "lambda", "\u03b1 err prob" = "alpha",
                "Power (1-\u03b2 err prob)" = "power", "\u03b2/\u03b1 ratio" = "error_ratio", "Numerator df" = "df1",
                "Denominator df" = "df2"),
  chi2_generic = c("Noncentrality parameter \u03bb" = "lambda", "\u03b1 err prob" = "alpha",
                   "Power (1-\u03b2 err prob)" = "power", "\u03b2/\u03b1 ratio" = "error_ratio", "Df" = "df")
)

g1_exact_labels <- c("Df", "Total sample size", "Sample size group 1", "Sample size group 2")

# MR-DISPLAY: closed rounding interval of the displayed digits plus the
# default allowance (fixed and scientific notation).
g1_display_ok <- function(ours, shown) {
  if (!is.numeric(ours) || length(ours) != 1L || !is.finite(ours)) return(FALSE)
  m <- regmatches(shown, regexec("^(-?[0-9]*)(\\.([0-9]*))?([eE]([+-]?[0-9]+))?$", shown))[[1]]
  dec <- nchar(m[4])
  ex <- if (nzchar(m[6])) as.numeric(m[6]) else 0
  x <- as.numeric(shown)
  half <- 0.5 * 10^(-dec) * 10^ex
  allowance <- 8 * .Machine$double.eps * max(1, abs(x))
  ours >= x - half - allowance && ours <= x + half + allowance
}

g1_spec <- function(cs) {
  map <- g1_input_map[[cs$method]]
  ins <- cs$inputs[setdiff(names(cs$inputs), c("Tail(s)", "drawer mode"))]
  inputs <- list()
  for (lab in names(ins)) {
    nm <- map[[lab]]
    v <- as.numeric(ins[[lab]])
    inputs[[nm]] <- if (nm %in% c("n_total", "n1", "n2")) as.integer(v) else v
  }
  options <- if ("Tail(s)" %in% names(cs$inputs)) list(tails = tolower(cs$inputs[["Tail(s)"]])) else list()
  analysis_spec(cs$method, cs$mode, inputs = inputs, options = options)
}

g1_shown_critical <- function(res) {
  crit <- res$quantities$critical
  ncp <- res$quantities$ncp
  if (res$method_id %in% c("f_generic", "chi2_generic")) return(crit)
  if (identical(res$mode, "criterion") && res$method_id != "z_generic") {
    return(if (length(crit) == 1L) (if (ncp < 0) -crit else crit) else max(crit))
  }
  if (length(crit) == 1L) crit else if (ncp < 0) min(crit) else max(crit)
}

g1_outputs <- function(res) {
  q <- res$quantities
  s <- res$solved
  a <- res$achieved
  gs <- res$group_sizes
  crit <- g1_shown_critical(res)
  out <- list("Noncentrality parameter \u03b4" = if (res$method_id == "t_generic") s$delta else q$ncp,
              "Noncentrality parameter \u03bc" = s$mu, "Noncentrality parameter \u03bb" = s$lambda,
              "Critical t" = crit, "Critical z" = crit, "Critical F" = crit, "Critical \u03c7\u00b2" = crit,
              "Df" = q$df, "Total sample size" = s$n_total %||% q$n_total,
              "Sample size group 1" = if (!is.null(gs)) unname(gs[["n1"]]),
              "Sample size group 2" = if (!is.null(gs)) unname(gs[["n2"]]),
              "Actual power" = a$power, "Power (1-\u03b2 err prob)" = a$power, "\u03b1 err prob" = s$alpha,
              "\u03b2 err prob" = a$beta, "Effect size d" = s$d, "Effect size dz" = s$dz)
  out
}

g1_drawer_value <- function(cs) {
  num <- function(lab) as.numeric(cs$inputs[[lab]])
  if (cs$method == "t_means_one_sample") {
    return(effect_size(cs$method, "from_means", mean_h0 = num("Mean H0"), mean_h1 = num("Mean H1"),
                       sd = num("SD \u03c3"))$magnitude)
  }
  if (cs$method == "t_means_matched_pairs") {
    if (grepl("From differences", cs$inputs[["drawer mode"]], fixed = TRUE)) {
      return(effect_size(cs$method, "from_differences", mean_diff = num("Mean of difference"),
                         sd_diff = num("SD of difference"))$value)
    }
    return(effect_size(cs$method, "from_group_parameters", mean_1 = num("Mean group 1"), mean_2 = num("Mean group 2"),
                       sd_1 = num("SD group 1"), sd_2 = num("SD group 2"),
                       rho = num("Correlation between groups"))$magnitude)
  }
  if (grepl("n1 != n2", cs$inputs[["drawer mode"]], fixed = TRUE)) {
    return(effect_size(cs$method, "common_sd", mean_1 = num("Mean group 1"), mean_2 = num("Mean group 2"),
                       sd = num("SD \u03c3 within each group"))$magnitude)
  }
  effect_size(cs$method, "two_sds", mean_1 = num("Mean group 1"), mean_2 = num("Mean group 2"),
              sd_1 = num("SD \u03c3 group 1"), sd_2 = num("SD \u03c3 group 2"))$magnitude
}

g1_check_case <- function(id) {
  cs <- g1_fixture_cases[[id]]
  if (identical(cs$mode, "effect_size_calculator")) {
    v <- g1_drawer_value(cs)
    for (lab in names(cs$outputs)) expect_true(g1_display_ok(v, cs$outputs[[lab]]), info = paste(id, lab, v))
    return(invisible())
  }
  if (!is.null(cs$message)) {
    # TD-0010: a refusal must be matched by a typed refusal without numbers.
    res <- tryCatch(calculate_power(g1_spec(cs)), inferguidR_error = function(e) e)
    code <- if (inherits(res, "inferguidR_error")) res$status$code else res$status$code
    expect_true(code %in% c("invalid_input", "infeasible"), info = paste(id, code))
    if (!inherits(res, "inferguidR_error")) {
      expect_identical(.ig_finite_paths(unclass(res)[c("solved", "achieved", "quantities")], ""), character(),
                       info = id)
    }
    return(invisible())
  }
  res <- calculate_power(g1_spec(cs))
  expect_identical(res$status$code, "ok", info = id)
  outs <- g1_outputs(res)
  for (lab in names(cs$outputs)) {
    v <- outs[[lab]]
    shown <- cs$outputs[[lab]]
    if (lab %in% g1_exact_labels) {
      expect_identical(as.double(v), as.numeric(shown), info = paste(id, lab))
    } else {
      expect_true(g1_display_ok(v, shown), info = sprintf("%s %s: G*Power %s, engine %.12g", id, lab, shown,
                                                          if (is.null(v)) NA_real_ else v))
    }
  }
}

g1_cases_of <- function(method) names(g1_fixture_cases)[vapply(g1_fixture_cases, `[[`, "", "method") == method]

for (m in unique(vapply(g1_fixture_cases, `[[`, "", "method"))) {
  local({
    method <- m
    test_that(sprintf("curated G*Power fixtures agree: %s", method), {
      ids <- g1_cases_of(method)
      for (id in ids) g1_check_case(id)
      modes <- unique(vapply(g1_fixture_cases[ids], `[[`, "", "mode"))
      offered <- .ig_get_method(method)$modes
      expect_true(all(offered %in% modes), info = paste(method, "every offered mode has a curated case"))
    })
  })
}

test_that("the embedded displayed strings equal the fixture records (development checkout)", {
  root <- normalizePath(testthat::test_path("..", ".."), mustWork = FALSE)
  fx <- file.path(root, "validation", "fixtures")
  skip_if_not(dir.exists(fx), "validation fixtures not available (installed package)")
  files <- c(list.files(file.path(fx, "emulator_capture"), pattern = "\\.json$", full.names = TRUE),
             list.files(file.path(fx, "emulator_capture_supplement"), pattern = "\\.json$", full.names = TRUE,
                        recursive = TRUE))
  index <- stats::setNames(files, sub("\\.json$", "", basename(files)))
  for (id in names(g1_fixture_cases)) {
    cs <- g1_fixture_cases[[id]]
    expect_true(id %in% names(index), info = id)
    rec <- jsonlite::parse_json(paste(readLines(index[[id]], encoding = "UTF-8", warn = FALSE), collapse = "\n"))
    outs <- Filter(function(o) !startsWith(o$name, "message"), rec$outputs)
    shown <- stats::setNames(vapply(outs, function(o) o$displayed, ""), vapply(outs, function(o) o$label_as_displayed, ""))
    shown <- shown[names(cs$outputs)]
    expect_identical(unname(shown), unname(cs$outputs), info = id)
    ins <- c(rec$inputs, rec$effect_size_drawer$inputs)
    vals <- stats::setNames(vapply(ins, function(i) i$displayed_in_reference %||% i$value_as_entered, ""),
                            vapply(ins, function(i) i$label_as_displayed, ""))
    main <- setdiff(names(cs$inputs), "drawer mode")
    expect_identical(unname(vals[main]), unname(cs$inputs[main]), info = id)
    expect_identical(rec$method$analysis_type, cs$mode, info = id)
    if (!is.null(cs$message)) {
      msg <- Filter(function(o) startsWith(o$name, "message"), rec$outputs)
      expect_identical(msg[[1]]$displayed, cs$message, info = id)
    }
  }
})
