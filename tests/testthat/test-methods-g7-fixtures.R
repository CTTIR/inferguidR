# Group G7: curated G*Power fixture comparisons (fast subset) -----------------------------
#
# A fast, offline subset of the comparisons of validation/compare/<method_id>/
# (full runs: validation/tools/compare_fixtures.R) for the four correlation z
# tests, at least one record per analysis type and option branch. Every case is
# a fixture record of validation/fixtures/emulator_capture/ or
# emulator_capture_supplement/g7/, cited by its record id; the displayed
# strings are copies of those records (G*Power 3.1.9.7 under Wine, evidence
# kind emulator_capture, decisions D22/D23). The last test checks the copies
# against the fixture files in a development checkout. Every case listed here
# agrees in the full comparison; the disagreeing records and their triage are
# in validation/compare/<method_id>/summary.md.
#
# Match rules (validation/tolerances.md, fixed before comparison): integers
# exact (MR-EXACT); displayed reals within the closed rounding interval of the
# last displayed digit plus 8 * eps * max(1, |x|) (MR-DISPLAY). G*Power
# refusals ('Input error') must be matched by a typed refusal without numbers
# (TD-0010). Display conventions (validation/compare/<method_id>/mapping.R):
# one critical z with the sign of the effect; tetrachoric two-tailed critical r
# truncated to [-1, 1] (exact mode) or the Bonett-Price interval under H0
# (approximation mode); the 'no solution' sentinel H1 corr = 2 of the dependent
# correlations is matched by the engine's status infeasible.

ref_reset_registry()

g7_fixture_cases <- list(
  "GPW-Z-CORR2-APR-B001-O01--r01" = list(method = "z_correlations_two_independent", mode = "a_priori",
    inputs = c("Tail(s)" = "One", "Effect size q" = "0.3", "\u03b1 err prob" = "0.05", "Power (1-\u03b2 err prob)" = "0.80", "Allocation ratio N2/N1" = "1"),
    outputs = c("Critical z" = "1.6448536", "Sample size group 1" = "141", "Sample size group 2" = "141", "Total sample size" = "282", "Actual power" = "0.8015397")),
  "GPW-Z-CORR2-APR-B001-D01--r01" = list(method = "z_correlations_two_independent", mode = "a_priori",
    inputs = c("Tail(s)" = "One", "Effect size q" = "0.01", "\u03b1 err prob" = "0.001", "Power (1-\u03b2 err prob)" = "0.99", "Allocation ratio N2/N1" = "0.2"),
    outputs = c("Critical z" = "3.0902323", "Sample size group 1" = "1760374", "Sample size group 2" = "352075", "Total sample size" = "2112449", "Actual power" = "0.9900001")),
  "GPW-Z-CORR2-APR-B002-E01--r01" = list(method = "z_correlations_two_independent", mode = "a_priori",
    inputs = c("Tail(s)" = "Two", "Effect size q" = "2.5", "\u03b1 err prob" = "0.05", "Power (1-\u03b2 err prob)" = "0.80", "Allocation ratio N2/N1" = "5"),
    outputs = c("Critical z" = "1.9599640", "Sample size group 1" = "5", "Sample size group 2" = "22", "Total sample size" = "27", "Actual power" = "0.9196916")),
  "SUP-G7-Z-CORR2-APR-A01--r01" = list(method = "z_correlations_two_independent", mode = "a_priori",
    inputs = c("Tail(s)" = "One", "Effect size q" = "0.5", "\u03b1 err prob" = "0.05", "Power (1-\u03b2 err prob)" = "0.80", "Allocation ratio N2/N1" = "1"),
    outputs = c("Critical z" = "1.6448536", "Sample size group 1" = "53", "Sample size group 2" = "53", "Total sample size" = "106", "Actual power" = "0.8037649")),
  "SUP-G7-Z-CORR2-APR-A03--r01" = list(method = "z_correlations_two_independent", mode = "a_priori",
    inputs = c("Tail(s)" = "One", "Effect size q" = "0.5", "\u03b1 err prob" = "0.05", "Power (1-\u03b2 err prob)" = "0.80", "Allocation ratio N2/N1" = "1.0001"),
    outputs = c("Critical z" = "1.6448536", "Sample size group 1" = "52", "Sample size group 2" = "53", "Total sample size" = "105", "Actual power" = "0.8002426")),
  "GPW-Z-CORR2-POS-B002-O01--r01" = list(method = "z_correlations_two_independent", mode = "post_hoc",
    inputs = c("Tail(s)" = "Two", "Effect size q" = "-0.4028126", "\u03b1 err prob" = "0.05", "Sample size group 1" = "260", "Sample size group 2" = "51"),
    outputs = c("Critical z" = "-1.9599640", "Power (1-\u03b2 err prob)" = "0.7263518")),
  "GPW-Z-CORR2-POS-B001-E01--r01" = list(method = "z_correlations_two_independent", mode = "post_hoc",
    inputs = c("Tail(s)" = "One", "Effect size q" = "-0.3", "\u03b1 err prob" = "0.05", "Sample size group 1" = "4", "Sample size group 2" = "4"),
    outputs = c("Critical z" = "-1.6448536", "Power (1-\u03b2 err prob)" = "0.0759687")),
  "SUP-G7-Z-CORR2-POS-A01--r01" = list(method = "z_correlations_two_independent", mode = "post_hoc",
    inputs = c("Tail(s)" = "One", "Effect size q" = "0.5", "\u03b1 err prob" = "0.05", "Sample size group 1" = "3", "Sample size group 2" = "50"),
    outputs = c("Critical z" = "1.6448536", "Power (1-\u03b2 err prob)" = "0.0500000")),
  "GPW-Z-CORR2-SEN-B001-O01--r01" = list(method = "z_correlations_two_independent", mode = "sensitivity",
    inputs = c("Tail(s)" = "One", "\u03b1 err prob" = "0.05", "Power (1-\u03b2 err prob)" = "0.80", "Sample size group 1" = "100", "Sample size group 2" = "100"),
    outputs = c("Critical z" = "1.6448536", "Effect size q" = "0.3570370")),
  "GPW-Z-CORR2-SEN-B002-D01--r01" = list(method = "z_correlations_two_independent", mode = "sensitivity",
    inputs = c("Tail(s)" = "Two", "\u03b1 err prob" = "0.001", "Power (1-\u03b2 err prob)" = "0.99", "Sample size group 1" = "100000", "Sample size group 2" = "100000"),
    outputs = c("Critical z" = "3.2905267", "Effect size q" = "0.0251198")),
  "GPW-Z-CORR2-CRI-B001-D01--r01" = list(method = "z_correlations_two_independent", mode = "criterion",
    inputs = c("Tail(s)" = "One", "Effect size q" = "0.01", "Power (1-\u03b2 err prob)" = "0.9999", "Sample size group 1" = "100000", "Sample size group 2" = "100000"),
    outputs = c("Critical z" = "-1.4829820", "\u03b1 err prob" = "0.9309604")),
  "GPW-Z-CORR2-CRI-B002-E01--r01" = list(method = "z_correlations_two_independent", mode = "criterion",
    inputs = c("Tail(s)" = "Two", "Effect size q" = "2.5", "Power (1-\u03b2 err prob)" = "0.95", "Sample size group 1" = "4", "Sample size group 2" = "4"),
    outputs = c("Critical z" = "0.2903912", "\u03b1 err prob" = "0.7715169")),
  "GPW-Z-CORR2-COM-B001-E01--r01" = list(method = "z_correlations_two_independent", mode = "compromise",
    inputs = c("Tail(s)" = "One", "Effect size q" = "2.5", "\u03b2/\u03b1 ratio" = "4", "Sample size group 1" = "4", "Sample size group 2" = "4"),
    outputs = c("Critical z" = "1.3659499", "\u03b1 err prob" = "0.0859773", "\u03b2 err prob" = "0.3439094", "Power (1-\u03b2 err prob)" = "0.6560906")),
  "GPW-Z-CORR2-COM-B002-O01--r01" = list(method = "z_correlations_two_independent", mode = "compromise",
    inputs = c("Tail(s)" = "Two", "Effect size q" = "0.3", "\u03b2/\u03b1 ratio" = "1", "Sample size group 1" = "100", "Sample size group 2" = "100"),
    outputs = c("Critical z" = "1.2671326", "\u03b1 err prob" = "0.2051079", "\u03b2 err prob" = "0.2051079", "Power (1-\u03b2 err prob)" = "0.7948921")),
  "GPW-Z-CORR2-APR-B001-I03--r01" = list(method = "z_correlations_two_independent", mode = "a_priori",
    inputs = c("Tail(s)" = "One", "Effect size q" = "0.3", "\u03b1 err prob" = "0.05", "Power (1-\u03b2 err prob)" = "0.04", "Allocation ratio N2/N1" = "1"),
    outputs = c("Critical z" = "?", "Sample size group 1" = "?", "Sample size group 2" = "?", "Total sample size" = "?", "Actual power" = "?"),
    message = "Power must be greater than or equal to alpha."),
  "GPW-Z-CORR2-COM-B001-I02--r01" = list(method = "z_correlations_two_independent", mode = "compromise",
    inputs = c("Tail(s)" = "One", "Effect size q" = "0.3", "\u03b2/\u03b1 ratio" = "1", "Sample size group 1" = "1", "Sample size group 2" = "1"),
    outputs = c("Critical z" = "?", "\u03b1 err prob" = "?", "\u03b2 err prob" = "?", "Power (1-\u03b2 err prob)" = "?"),
    message = "Total sample size is outside the admissible range (3, 2147483646)"),
  "GPW-Z-DEPCI-APR-B001-O01--r01" = list(method = "z_correlations_two_dependent_common_index", mode = "a_priori",
    inputs = c("Tail(s)" = "One", "H1 corr \u03c1_ac" = "0.2", "\u03b1 err prob" = "0.05", "Power (1-\u03b2 err prob)" = "0.8", "H0 corr \u03c1_ab" = "0.4", "Corr \u03c1_bc" = "0.5"),
    outputs = c("Critical z" = "-1.6448536", "Sample size" = "144", "Actual power" = "0.8011611")),
  "GPW-Z-DEPCI-APR-B004-E01--r01" = list(method = "z_correlations_two_dependent_common_index", mode = "a_priori",
    inputs = c("Tail(s)" = "Two", "H1 corr \u03c1_ac" = "-0.6", "\u03b1 err prob" = "0.05", "Power (1-\u03b2 err prob)" = "0.80", "H0 corr \u03c1_ab" = "0.6", "Corr \u03c1_bc" = "0"),
    outputs = c("Critical z" = "-1.9599640", "Sample size" = "12", "Actual power" = "0.8218221")),
  "GPW-Z-DEPCI-POS-B004-E01--r01" = list(method = "z_correlations_two_dependent_common_index", mode = "post_hoc",
    inputs = c("Tail(s)" = "Two", "H1 corr \u03c1_ac" = "0.7", "\u03b1 err prob" = "0.05", "Sample size" = "4", "H0 corr \u03c1_ab" = "0.6", "Corr \u03c1_bc" = "0"),
    outputs = c("Critical z" = "1.9599640", "Power (1-\u03b2 err prob)" = "0.0488391")),
  "GPW-Z-DEPCI-SEN-B002-O01--r01" = list(method = "z_correlations_two_dependent_common_index", mode = "sensitivity",
    inputs = c("Tail(s)" = "One", "Effect direction" = "\u03c1_ac \u2264 \u03c1_ab", "\u03b1 err prob" = "0.05", "Power (1-\u03b2 err prob)" = "0.8", "Sample size" = "144", "H0 corr \u03c1_ab" = "0.4", "Corr \u03c1_bc" = "-0.6"),
    outputs = c("Critical z" = "-1.6448536", "H1 corr \u03c1_ac" = "0.0477022")),
  "GPW-Z-DEPCI-SEN-B006-D01--r01" = list(method = "z_correlations_two_dependent_common_index", mode = "sensitivity",
    inputs = c("Tail(s)" = "Two", "Effect direction" = "\u03c1_ac \u2265 \u03c1_ab", "\u03b1 err prob" = "0.001", "Power (1-\u03b2 err prob)" = "0.99", "Sample size" = "100000", "H0 corr \u03c1_ab" = "0.4", "Corr \u03c1_bc" = "0.9"),
    outputs = c("Critical z" = "3.2905267", "H1 corr \u03c1_ac" = "0.4072253")),
  "GPW-Z-DEPCI-SEN-B003-E01--r01" = list(method = "z_correlations_two_dependent_common_index", mode = "sensitivity",
    inputs = c("Tail(s)" = "One", "Effect direction" = "\u03c1_ac \u2265 \u03c1_ab", "\u03b1 err prob" = "0.05", "Power (1-\u03b2 err prob)" = "0.8", "Sample size" = "144", "H0 corr \u03c1_ab" = "0.4", "Corr \u03c1_bc" = "-0.6"),
    outputs = c("Critical z" = "1.6448536", "H1 corr \u03c1_ac" = "2.0000000"),
    message = "No solution found for the requested effect size!  This is indicated by a (nonsensical) correlation value of 2. The most probable cause is that there is no value leading to a valid (i.e. positive semi-definite) correlation matrix."),
  "SUP-G7-Z-DEPCI-SEN-A01--r01" = list(method = "z_correlations_two_dependent_common_index", mode = "sensitivity",
    inputs = c("Tail(s)" = "Two", "Effect direction" = "\u03c1_ac \u2265 \u03c1_ab", "\u03b1 err prob" = "0.05", "Power (1-\u03b2 err prob)" = "0.80", "Sample size" = "5", "H0 corr \u03c1_ab" = "0.4", "Corr \u03c1_bc" = "0.5"),
    outputs = c("Critical z" = "1.9599640", "H1 corr \u03c1_ac" = "0.9877916")),
  "GPW-Z-DEPCI-CRI-B001-E01--r01" = list(method = "z_correlations_two_dependent_common_index", mode = "criterion",
    inputs = c("Tail(s)" = "One", "H1 corr \u03c1_ac" = "-0.6", "Power (1-\u03b2 err prob)" = "0.95", "Sample size" = "4", "H0 corr \u03c1_ab" = "0.6", "Corr \u03c1_bc" = "0"),
    outputs = c("Critical z" = "0.5285054", "\u03b1 err prob" = "0.7014257")),
  "GPW-Z-DEPCI-CRI-B004-O01--r01" = list(method = "z_correlations_two_dependent_common_index", mode = "criterion",
    inputs = c("Tail(s)" = "Two", "H1 corr \u03c1_ac" = "0.2", "Power (1-\u03b2 err prob)" = "0.80", "Sample size" = "144", "H0 corr \u03c1_ab" = "0.4", "Corr \u03c1_bc" = "0.5"),
    outputs = c("Critical z" = "-1.6489129", "\u03b1 err prob" = "0.0991655")),
  "GPW-Z-DEPCI-COM-B001-D01--r01" = list(method = "z_correlations_two_dependent_common_index", mode = "compromise",
    inputs = c("Tail(s)" = "One", "H1 corr \u03c1_ac" = "0.41", "\u03b2/\u03b1 ratio" = "0.01", "Sample size" = "100000", "H0 corr \u03c1_ab" = "0.4", "Corr \u03c1_bc" = "0.9"),
    outputs = c("Critical z" = "3.3312006", "\u03b1 err prob" = "0.0004323614", "\u03b2 err prob" = "4.323614e-06", "Power (1-\u03b2 err prob)" = "0.9999957")),
  "GPW-Z-DEPCI-APR-B001-I02--r01" = list(method = "z_correlations_two_dependent_common_index", mode = "a_priori",
    inputs = c("Tail(s)" = "One", "H1 corr \u03c1_ac" = "0.4", "\u03b1 err prob" = "0.05", "Power (1-\u03b2 err prob)" = "0.80", "H0 corr \u03c1_ab" = "0.4", "Corr \u03c1_bc" = "0.5"),
    outputs = c("Critical z" = "?", "Sample size" = "?", "Actual power" = "?"),
    message = "Effect size must not be zero, i.e. check that |r_cd - r_ab| > 1e-6."),
  "GPW-Z-DEPNCI-APR-B001-O01--r01" = list(method = "z_correlations_two_dependent_no_common_index", mode = "a_priori",
    inputs = c("Tail(s)" = "One", "H1 corr \u03c1_cd" = "0.2", "\u03b1 err prob" = "0.05", "Power (1-\u03b2 err prob)" = "0.8", "H0 corr \u03c1_ab" = "0.1", "Corr \u03c1_ac" = "0.5", "Corr \u03c1_ad" = "0.4", "Corr \u03c1_bc" = "-0.4", "Corr \u03c1_bd" = "0.8"),
    outputs = c("Critical z" = "1.6448536", "Sample size" = "886", "Actual power" = "0.8001562")),
  "GPW-Z-DEPNCI-APR-B004-D01--r01" = list(method = "z_correlations_two_dependent_no_common_index", mode = "a_priori",
    inputs = c("Tail(s)" = "Two", "H1 corr \u03c1_cd" = "0.31", "\u03b1 err prob" = "0.001", "Power (1-\u03b2 err prob)" = "0.99", "H0 corr \u03c1_ab" = "0.3", "Corr \u03c1_ac" = "0.2", "Corr \u03c1_ad" = "0.2", "Corr \u03c1_bc" = "0.2", "Corr \u03c1_bd" = "0.2"),
    outputs = c("Critical z" = "3.2905267", "Sample size" = "494549", "Actual power" = "0.9900001")),
  "GPW-Z-DEPNCI-POS-B001-E01--r01" = list(method = "z_correlations_two_dependent_no_common_index", mode = "post_hoc",
    inputs = c("Tail(s)" = "One", "H1 corr \u03c1_cd" = "-0.5", "\u03b1 err prob" = "0.05", "Sample size" = "4", "H0 corr \u03c1_ab" = "0.1", "Corr \u03c1_ac" = "0", "Corr \u03c1_ad" = "0", "Corr \u03c1_bc" = "0", "Corr \u03c1_bd" = "0"),
    outputs = c("Critical z" = "-1.6448536", "Power (1-\u03b2 err prob)" = "0.1179123")),
  "GPW-Z-DEPNCI-SEN-B005-E01--r01" = list(method = "z_correlations_two_dependent_no_common_index", mode = "sensitivity",
    inputs = c("Tail(s)" = "Two", "Effect direction" = "\u03c1_cd \u2264 \u03c1_ab", "\u03b1 err prob" = "0.05", "Power (1-\u03b2 err prob)" = "0.95", "Sample size" = "4", "H0 corr \u03c1_ab" = "0.1", "Corr \u03c1_ac" = "0", "Corr \u03c1_ad" = "0", "Corr \u03c1_bc" = "0", "Corr \u03c1_bd" = "0"),
    outputs = c("Critical z" = "-1.9599640", "H1 corr \u03c1_cd" = "-0.9999088")),
  "GPW-Z-DEPNCI-SEN-B003-O01--r01" = list(method = "z_correlations_two_dependent_no_common_index", mode = "sensitivity",
    inputs = c("Tail(s)" = "One", "Effect direction" = "\u03c1_cd \u2265 \u03c1_ab", "\u03b1 err prob" = "0.05", "Power (1-\u03b2 err prob)" = "0.80", "Sample size" = "500", "H0 corr \u03c1_ab" = "0.1", "Corr \u03c1_ac" = "0.5", "Corr \u03c1_ad" = "0.4", "Corr \u03c1_bc" = "-0.4", "Corr \u03c1_bd" = "0.8"),
    outputs = c("Critical z" = "1.6448536", "H1 corr \u03c1_cd" = "2.0000000"),
    message = "No solution found for the requested effect size!  This is indicated by a (nonsensical) correlation value of 2. The most probable cause is that there is no value leading to a valid (i.e. positive semi-definite) correlation matrix."),
  "SUP-G7-Z-DEPNCI-SEN-A02--r01" = list(method = "z_correlations_two_dependent_no_common_index", mode = "sensitivity",
    inputs = c("Tail(s)" = "One", "Effect direction" = "\u03c1_cd \u2265 \u03c1_ab", "\u03b1 err prob" = "0.05", "Power (1-\u03b2 err prob)" = "0.80", "Sample size" = "540", "H0 corr \u03c1_ab" = "0.1", "Corr \u03c1_ac" = "0.5", "Corr \u03c1_ad" = "0.4", "Corr \u03c1_bc" = "-0.4", "Corr \u03c1_bd" = "0.8"),
    outputs = c("Critical z" = "1.6448536", "H1 corr \u03c1_cd" = "0.2275172")),
  "GPW-Z-DEPNCI-CRI-B001-D01--r01" = list(method = "z_correlations_two_dependent_no_common_index", mode = "criterion",
    inputs = c("Tail(s)" = "One", "H1 corr \u03c1_cd" = "0.31", "Power (1-\u03b2 err prob)" = "0.9999", "Sample size" = "100000", "H0 corr \u03c1_ab" = "0.3", "Corr \u03c1_ac" = "0.2", "Corr \u03c1_ad" = "0.2", "Corr \u03c1_bc" = "0.2", "Corr \u03c1_bd" = "0.2"),
    outputs = c("Critical z" = "-1.1938039", "\u03b1 err prob" = "0.8837227")),
  "GPW-Z-DEPNCI-COM-B004-O01--r01" = list(method = "z_correlations_two_dependent_no_common_index", mode = "compromise",
    inputs = c("Tail(s)" = "Two", "H1 corr \u03c1_cd" = "0.2", "\u03b2/\u03b1 ratio" = "1", "Sample size" = "500", "H0 corr \u03c1_ab" = "0.1", "Corr \u03c1_ac" = "0.5", "Corr \u03c1_ad" = "0.4", "Corr \u03c1_bc" = "-0.4", "Corr \u03c1_bd" = "0.8"),
    outputs = c("Critical z" = "1.1699355", "\u03b1 err prob" = "0.2420269", "\u03b2 err prob" = "0.2420269", "Power (1-\u03b2 err prob)" = "0.7579731")),
  "GPW-Z-DEPNCI-POS-B001-I02--r01" = list(method = "z_correlations_two_dependent_no_common_index", mode = "post_hoc",
    inputs = c("Tail(s)" = "One", "H1 corr \u03c1_cd" = "0.2", "\u03b1 err prob" = "0.05", "Sample size" = "500", "H0 corr \u03c1_ab" = "0.1", "Corr \u03c1_ac" = "0.5", "Corr \u03c1_ad" = "0.4", "Corr \u03c1_bc" = "-0.4", "Corr \u03c1_bd" = "0.9"),
    outputs = c("Critical z" = "?", "Power (1-\u03b2 err prob)" = "?"),
    message = "Under H0: The correlation matrix is not valid, that is, not positive semi-definite."),
  "GPW-Z-TETRA-APR-B001-O01--r01" = list(method = "z_tetrachoric_correlation", mode = "a_priori", computation = "Use exact computation (Brown and Benedetti, 1977)",
    inputs = c("Tail(s)" = "One", "H1 corr \u03c1" = "0.2399846", "\u03b1 err prob" = "0.05", "Power (1-\u03b2 err prob)" = "0.95", "H0 corr \u03c1" = "0", "Marginal prop x" = "0.6019313", "Marginal prop y" = "0.5815451"),
    outputs = c("Critical z" = "1.6448536", "Total sample size" = "463", "Actual power" = "0.9503705", "H1 corr \u03c1" = "0.2399846", "H0 corr \u03c1" = "-1.110223e-16", "Critical r lwr" = "0.1224841", "Critical r upr" = "0.1224841", "Std err r" = "0.0744650")),
  "GPW-Z-TETRA-APR-B010-D01--r01" = list(method = "z_tetrachoric_correlation", mode = "a_priori", computation = "Use approximation (Bonett and Price, 2005)",
    inputs = c("Tail(s)" = "Two", "H1 corr \u03c1" = "0.22", "\u03b1 err prob" = "0.001", "Power (1-\u03b2 err prob)" = "0.99", "H0 corr \u03c1" = "0.2", "Marginal prop x" = "0.5", "Marginal prop y" = "0.5"),
    outputs = c("Critical z" = "3.2905267", "Total sample size" = "182821", "Actual power" = "0.9900002", "H1 corr \u03c1*" = "0.2200000", "H0 corr \u03c1*" = "0.2000000", "Critical r* lwr" = "0.1882279", "Critical r* upr" = "0.2117200", "Std err r*" = "0.003569806")),
  "SUP-G7-Z-TETRA-APR-A02--r01" = list(method = "z_tetrachoric_correlation", mode = "a_priori", computation = "Use exact computation (Brown and Benedetti, 1977)",
    inputs = c("Tail(s)" = "Two", "H1 corr \u03c1" = "0.35", "\u03b1 err prob" = "0.05", "Power (1-\u03b2 err prob)" = "0.8", "H0 corr \u03c1" = "0.05", "Marginal prop x" = "0.4", "Marginal prop y" = "0.25"),
    outputs = c("Critical z" = "1.9599640", "Total sample size" = "245", "Actual power" = "0.8006253", "H1 corr \u03c1" = "0.3500000", "H0 corr \u03c1" = "0.0500000", "Critical r lwr" = "-0.1656324", "Critical r upr" = "0.2656324", "Std err r" = "0.1100186")),
  "GPW-Z-TETRA-POS-B007-E01--r01" = list(method = "z_tetrachoric_correlation", mode = "post_hoc", computation = "Use exact computation (Brown and Benedetti, 1977)",
    inputs = c("Tail(s)" = "Two", "H1 corr \u03c1" = "-0.3", "\u03b1 err prob" = "0.05", "Total sample size" = "10", "H0 corr \u03c1" = "0", "Marginal prop x" = "0.5", "Marginal prop y" = "0.5"),
    outputs = c("Critical z" = "-1.9599640", "Power (1-\u03b2 err prob)" = "0.0767417", "H1 corr \u03c1" = "-0.3000000", "H0 corr \u03c1" = "0", "Critical r lwr" = "-0.9735718", "Critical r upr" = "0.9735718", "Std err r" = "0.4967294")),
  "SUP-G7-Z-TETRA-POS-A01--r01" = list(method = "z_tetrachoric_correlation", mode = "post_hoc", computation = "Use exact computation (Brown and Benedetti, 1977)",
    inputs = c("Tail(s)" = "Two", "H1 corr \u03c1" = "0.7", "\u03b1 err prob" = "0.05", "Total sample size" = "8", "H0 corr \u03c1" = "0.5", "Marginal prop x" = "0.3", "Marginal prop y" = "0.6"),
    outputs = c("Critical z" = "1.9599640", "Power (1-\u03b2 err prob)" = "0.0259036", "H1 corr \u03c1" = "0.7000000", "H0 corr \u03c1" = "0.5000000", "Critical r lwr" = "-0.4719314", "Critical r upr" = "1.0000000", "Std err r" = "0.4958925")),
  "SUP-G7-Z-TETRA-POS-A02--r01" = list(method = "z_tetrachoric_correlation", mode = "post_hoc", computation = "Use approximation (Bonett and Price, 2005)",
    inputs = c("Tail(s)" = "Two", "H1 corr \u03c1" = "0.4", "\u03b1 err prob" = "0.05", "Total sample size" = "150", "H0 corr \u03c1" = "0.1", "Marginal prop x" = "0.3", "Marginal prop y" = "0.7"),
    outputs = c("Critical z" = "1.9599640", "Power (1-\u03b2 err prob)" = "0.5986893", "H1 corr \u03c1*" = "0.3983425", "H0 corr \u03c1*" = "0.0956036", "Critical r* lwr" = "-0.1729758", "Critical r* upr" = "0.3520909", "Std err r*" = "0.1368367")),
  "SUP-G7-Z-TETRA-POS-A05--r01" = list(method = "z_tetrachoric_correlation", mode = "post_hoc", computation = "Use exact computation (Brown and Benedetti, 1977)",
    inputs = c("Tail(s)" = "One", "H1 corr \u03c1" = "0.5", "\u03b1 err prob" = "0.05", "Total sample size" = "100", "H0 corr \u03c1" = "0", "Marginal prop x" = "0.05", "Marginal prop y" = "0.95"),
    outputs = c("Critical z" = "1.6448536", "Power (1-\u03b2 err prob)" = "0.4016530", "H1 corr \u03c1" = "0.5000000", "H0 corr \u03c1" = "6.767074e-16", "Critical r lwr" = "0.7345195", "Critical r upr" = "0.7345195", "Std err r" = "0.4465561")),
  "GPW-Z-TETRA-SEN-B009-O01--r01" = list(method = "z_tetrachoric_correlation", mode = "sensitivity", computation = "Use exact computation (Brown and Benedetti, 1977)",
    inputs = c("Tail(s)" = "Two", "Effect direction" = "\u03c1 H1 \u2265 \u03c1 H0", "\u03b1 err prob" = "0.05", "Power (1-\u03b2 err prob)" = "0.80", "Total sample size" = "200", "H0 corr \u03c1" = "0", "Marginal prop x" = "0.5", "Marginal prop y" = "0.5"),
    outputs = c("Critical z" = "1.9599640", "H1 corr \u03c1" = "0.3049745", "H0 corr \u03c1" = "0", "Critical r lwr" = "-0.2176973", "Critical r upr" = "0.2176973", "Std err r" = "0.1110721")),
  "SUP-G7-Z-TETRA-SEN-A02--r01" = list(method = "z_tetrachoric_correlation", mode = "sensitivity", computation = "Use exact computation (Brown and Benedetti, 1977)",
    inputs = c("Tail(s)" = "Two", "Effect direction" = "\u03c1 H1 \u2264 \u03c1 H0", "\u03b1 err prob" = "0.05", "Power (1-\u03b2 err prob)" = "0.8", "Total sample size" = "150", "H0 corr \u03c1" = "0.1", "Marginal prop x" = "0.3", "Marginal prop y" = "0.7"),
    outputs = c("Critical z" = "-1.9599640", "H1 corr \u03c1" = "-0.2886073", "H0 corr \u03c1" = "0.1000000", "Critical r lwr" = "-0.1777392", "Critical r upr" = "0.3777392", "Std err r" = "0.1417063")),
  "GPW-Z-TETRA-CRI-B007-O01--r01" = list(method = "z_tetrachoric_correlation", mode = "criterion", computation = "Use exact computation (Brown and Benedetti, 1977)",
    inputs = c("Tail(s)" = "Two", "H1 corr \u03c1" = "0.3", "Power (1-\u03b2 err prob)" = "0.80", "Total sample size" = "200", "H0 corr \u03c1" = "0", "Marginal prop x" = "0.5", "Marginal prop y" = "0.5"),
    outputs = c("Critical z" = "1.9133436", "\u03b1 err prob" = "0.0557041", "H1 corr \u03c1" = "0.3000000", "H0 corr \u03c1" = "0", "Critical r lwr" = "-0.2125190", "Critical r upr" = "0.2125190", "Std err r" = "0.1110721")),
  "SUP-G7-Z-TETRA-CRI-A02--r01" = list(method = "z_tetrachoric_correlation", mode = "criterion", computation = "Use approximation (Bonett and Price, 2005)",
    inputs = c("Tail(s)" = "Two", "H1 corr \u03c1" = "0.4", "Power (1-\u03b2 err prob)" = "0.8", "Total sample size" = "150", "H0 corr \u03c1" = "0.1", "Marginal prop x" = "0.3", "Marginal prop y" = "0.7"),
    outputs = c("Critical z" = "1.3629677", "\u03b1 err prob" = "0.1728927", "H1 corr \u03c1*" = "0.3983425", "H0 corr \u03c1*" = "0.0956036", "Critical r* lwr" = "-0.0919565", "Critical r* upr" = "0.2771987", "Std err r*" = "0.1368367")),
  "GPW-Z-TETRA-COM-B004-D01--r01" = list(method = "z_tetrachoric_correlation", mode = "compromise", computation = "Use approximation (Bonett and Price, 2005)",
    inputs = c("Tail(s)" = "One", "H1 corr \u03c1" = "0.22", "\u03b2/\u03b1 ratio" = "0.01", "Total sample size" = "100000", "H0 corr \u03c1" = "0.2", "Marginal prop x" = "0.5", "Marginal prop y" = "0.5"),
    outputs = c("Critical z" = "1.1426459", "\u03b1 err prob" = "0.1265928", "\u03b2 err prob" = "0.001265928", "Power (1-\u03b2 err prob)" = "0.9987341", "H1 corr \u03c1*" = "0.2200000", "H0 corr \u03c1*" = "0.2000000", "Critical r* lwr" = "0.2055153", "Critical r* upr" = "0.2055153", "Std err r*" = "0.004826781")),
  "SUP-G7-Z-TETRA-COM-A01--r01" = list(method = "z_tetrachoric_correlation", mode = "compromise", computation = "Use exact computation (Brown and Benedetti, 1977)",
    inputs = c("Tail(s)" = "Two", "H1 corr \u03c1" = "0.4", "\u03b2/\u03b1 ratio" = "2", "Total sample size" = "150", "H0 corr \u03c1" = "0.1", "Marginal prop x" = "0.3", "Marginal prop y" = "0.7"),
    outputs = c("Critical z" = "1.5190957", "\u03b1 err prob" = "0.1287384", "\u03b2 err prob" = "0.2574768", "Power (1-\u03b2 err prob)" = "0.7425232", "H1 corr \u03c1" = "0.4000000", "H0 corr \u03c1" = "0.1000000", "Critical r lwr" = "-0.1152654", "Critical r upr" = "0.3152654", "Std err r" = "0.1417063")),
  "GPW-Z-TETRA-APR-B001-I02--r01" = list(method = "z_tetrachoric_correlation", mode = "a_priori", computation = "Use exact computation (Brown and Benedetti, 1977)",
    inputs = c("Tail(s)" = "One", "H1 corr \u03c1" = "0", "\u03b1 err prob" = "0.05", "Power (1-\u03b2 err prob)" = "0.80", "H0 corr \u03c1" = "0", "Marginal prop x" = "0.5", "Marginal prop y" = "0.5"),
    outputs = c("Critical z" = "?", "Total sample size" = "?", "Actual power" = "?", "H1 corr \u03c1" = "?", "H0 corr \u03c1" = "?", "Critical r lwr" = "?", "Critical r upr" = "?", "Std err r" = "?"),
    message = "Effect size to small (|H0 rho - H1 rho| < 1e-6)")
)

g7_label_map <- list(
  z_correlations_two_independent = c("Effect size q" = "q", "\u03b1 err prob" = "alpha", "Power (1-\u03b2 err prob)" = "power",
                                     "\u03b2/\u03b1 ratio" = "error_ratio", "Sample size group 1" = "n1",
                                     "Sample size group 2" = "n2", "Allocation ratio N2/N1" = "allocation_ratio"),
  z_correlations_two_dependent_common_index = c("H1 corr \u03c1_ac" = "rho_ac", "\u03b1 err prob" = "alpha",
                                                "Power (1-\u03b2 err prob)" = "power", "\u03b2/\u03b1 ratio" = "error_ratio",
                                                "Sample size" = "n_total", "H0 corr \u03c1_ab" = "rho_ab", "Corr \u03c1_bc" = "rho_bc"),
  z_correlations_two_dependent_no_common_index = c("H1 corr \u03c1_cd" = "rho_cd", "\u03b1 err prob" = "alpha",
                                                   "Power (1-\u03b2 err prob)" = "power", "\u03b2/\u03b1 ratio" = "error_ratio",
                                                   "Sample size" = "n_total", "H0 corr \u03c1_ab" = "rho_ab",
                                                   "Corr \u03c1_ac" = "rho_ac", "Corr \u03c1_ad" = "rho_ad",
                                                   "Corr \u03c1_bc" = "rho_bc", "Corr \u03c1_bd" = "rho_bd"),
  z_tetrachoric_correlation = c("H1 corr \u03c1" = "rho_h1", "\u03b1 err prob" = "alpha", "Power (1-\u03b2 err prob)" = "power",
                                "\u03b2/\u03b1 ratio" = "error_ratio", "Total sample size" = "n_total",
                                "H0 corr \u03c1" = "rho_h0", "Marginal prop x" = "marginal_x", "Marginal prop y" = "marginal_y")
)

g7_direction_codes <- c("\u03c1_ac \u2264 \u03c1_ab" = "rho_ac_le_rho_ab", "\u03c1_ac \u2265 \u03c1_ab" = "rho_ac_ge_rho_ab",
                        "\u03c1_cd \u2264 \u03c1_ab" = "rho_cd_le_rho_ab", "\u03c1_cd \u2265 \u03c1_ab" = "rho_cd_ge_rho_ab",
                        "\u03c1 H1 \u2264 \u03c1 H0" = "rho_h1_le_h0", "\u03c1 H1 \u2265 \u03c1 H0" = "rho_h1_ge_h0")

g7_case_spec <- function(cs) {
  map <- g7_label_map[[cs$method]]
  inputs <- list()
  options <- list()
  for (lab in names(cs$inputs)) {
    s <- cs$inputs[[lab]]
    if (lab == "Tail(s)") {
      options$tails <- tolower(s)
    } else if (lab == "Effect direction") {
      options$effect_direction <- g7_direction_codes[[s]]
    } else {
      nm <- map[[lab]]
      v <- as.numeric(s)
      inputs[[nm]] <- if (nm %in% c("n1", "n2", "n_total")) as.integer(v) else v
    }
  }
  if (!is.null(cs$computation)) {
    options$computation_of_correlation <- if (startsWith(cs$computation, "Use approximation")) "approx_bonett_price" else
      "exact_brown_benedetti"
  }
  list(inputs = inputs, options = options)
}

# Engine values by G*Power output label (the mapping decisions of
# validation/compare/<method_id>/mapping.R).
g7_engine_outputs <- function(cs, res) {
  q <- res$quantities
  mo <- res$method_outputs
  crit <- q$critical
  crit_z <- if (length(crit) == 1L) crit else if (isTRUE(q$ncp < 0)) min(crit) else max(crit)
  out <- list("Critical z" = crit_z, "Actual power" = res$achieved$power, "Power (1-\u03b2 err prob)" = res$achieved$power,
              "\u03b1 err prob" = res$solved$alpha, "\u03b2 err prob" = res$achieved$beta)
  if (cs$method == "z_correlations_two_independent") {
    out[["Sample size group 1"]] <- unname(res$group_sizes[["n1"]] %||% NA)
    out[["Sample size group 2"]] <- unname(res$group_sizes[["n2"]] %||% NA)
    out[["Total sample size"]] <- res$solved$n_total
    out[["Effect size q"]] <- res$solved$q
  } else if (cs$method == "z_tetrachoric_correlation") {
    approx <- identical(mo$computation, "approx_bonett_price")
    st <- if (approx) "*" else ""
    cr <- mo$critical_r
    lim <- if (length(cr) == 1L) c(cr, cr) else if (approx) mo$bonett_price_interval_h0 else c(max(-1, cr[1]), min(1, cr[2]))
    out[["Total sample size"]] <- res$solved$n_total
    out[[paste0("H1 corr \u03c1", st)]] <- mo$rho_h1_used
    out[[paste0("H0 corr \u03c1", st)]] <- mo$rho_h0_used
    out[[paste0("Critical r", st, " lwr")]] <- lim[1]
    out[[paste0("Critical r", st, " upr")]] <- lim[2]
    out[[paste0("Std err r", st)]] <- mo$se_h0
  } else {
    h1 <- if (cs$method == "z_correlations_two_dependent_common_index") "rho_ac" else "rho_cd"
    out[["Sample size"]] <- res$solved$n_total
    out[[paste0("H1 corr \u03c1_", substr(h1, 5, 6))]] <- res$solved[[h1]]
  }
  out
}

g7_display_ok <- function(ours, shown) {
  if (is.null(ours) || length(ours) != 1L || !is.finite(ours)) return(FALSE)
  x <- as.numeric(shown)
  if (grepl("e", shown, fixed = TRUE)) {
    mant <- sub("e.*$", "", shown)
    dec <- if (grepl(".", mant, fixed = TRUE)) nchar(sub("^[^.]*\\.", "", mant)) else 0
    half <- 0.5 * 10^(-dec) * 10^as.numeric(sub("^.*e", "", shown))
  } else {
    dec <- if (grepl(".", shown, fixed = TRUE)) nchar(sub("^[^.]*\\.", "", shown)) else 0
    half <- if (dec == 0 && !grepl(".", shown, fixed = TRUE) && x == round(x) && abs(x) > 1) 0 else 0.5 * 10^-dec
  }
  allowance <- 8 * .Machine$double.eps * max(1, abs(x))
  ours >= x - half - allowance && ours <= x + half + allowance
}

test_that("curated G7 fixture records agree with the engine (all modes and option branches)", {
  for (id in names(g7_fixture_cases)) {
    cs <- g7_fixture_cases[[id]]
    sp <- g7_case_spec(cs)
    if (!is.null(cs$message)) {
      code <- tryCatch({
        r <- calculate_power(analysis_spec(cs$method, cs$mode, inputs = sp$inputs, options = sp$options))
        r$status$code
      }, inferguidR_error = function(e) e$status$code)
      expect_true(code %in% c("invalid_input", "infeasible"), info = sprintf("%s (%s): engine %s", id, cs$message, code))
      next
    }
    res <- calculate_power(analysis_spec(cs$method, cs$mode, inputs = sp$inputs, options = sp$options))
    sentinel <- any(cs$outputs == "2.0000000") && cs$mode == "sensitivity" && cs$method != "z_tetrachoric_correlation"
    if (sentinel) {
      expect_identical(res$status$code, "infeasible", info = id)
      next
    }
    expect_identical(res$status$code, "ok", info = id)
    eng <- g7_engine_outputs(cs, res)
    for (lab in names(cs$outputs)) {
      expect_true(g7_display_ok(eng[[lab]], cs$outputs[[lab]]),
                  info = sprintf("%s %s: G*Power %s, engine %s", id, lab, cs$outputs[[lab]], format(eng[[lab]], digits = 12)))
    }
  }
})

test_that("the curated copies match the fixture files (development checkout only)", {
  root <- th_pkg_root()
  skip_if(is.null(root), "not a development checkout")
  fx <- file.path(root, "validation", "fixtures")
  skip_if_not(dir.exists(fx), "fixture records are not part of the built package")
  files <- c(list.files(file.path(fx, "emulator_capture"), pattern = "\\.json$", full.names = TRUE),
             list.files(file.path(fx, "emulator_capture_supplement", "g7"), pattern = "\\.json$", full.names = TRUE))
  idx <- stats::setNames(files, sub("\\.json$", "", basename(files)))
  for (id in names(g7_fixture_cases)) {
    cs <- g7_fixture_cases[[id]]
    expect_true(id %in% names(idx), info = id)
    rec <- jsonlite::parse_json(paste(readLines(idx[[id]], encoding = "UTF-8", warn = FALSE), collapse = "\n"))
    expect_identical(rec$method$method_id, cs$method, info = id)
    shown_in <- stats::setNames(vapply(rec$inputs, function(i) i$displayed_in_reference, ""),
                                vapply(rec$inputs, function(i) i$label_as_displayed, ""))
    expect_identical(shown_in[names(cs$inputs)], cs$inputs, info = id)
    vals <- Filter(function(o) !startsWith(o$name, "message"), rec$outputs)
    shown_out <- stats::setNames(vapply(vals, function(o) o$displayed, ""), vapply(vals, function(o) o$label_as_displayed, ""))
    expect_identical(shown_out[names(cs$outputs)], cs$outputs, info = id)
  }
})
