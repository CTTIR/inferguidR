# Protocol export (R/protocol-export.R). Group G12; method notes
# specs/method-notes/utilities.md, section "Protocol export".
#
# Expected texts are G*Power 3.1.9.7 protocol exports (emulator_capture under
# Wine) of the cited fixture records; the raw files are
# validation/raw/emulator/cases/<method>/<record_id>__protocol__01.txt with the
# SHA-256 recorded in each fixture. The corpus-wide comparison of every method
# is validation/compare/utilities/util_protocol/protocol-summary.md.

t2 <- "t_means_two_independent"
t2_res <- function(mode, tails = "two", ...) {
  calculate_power(analysis_spec(t2, mode, inputs = list(...), options = list(tails = tails)))
}
crlf <- function(...) paste0(paste(c(...), collapse = "\r\n"), "\r\n")
strip_end <- function(s) sub("[\r\n]+$", "", s)

test_that("a post hoc record is identical to the G*Power protocol export", {
  # GPW-UT-PROT-UTL-U002-UO01--r01, raw protocol SHA-256 9babde83...cc9042a3.
  res <- t2_res("post_hoc", d = 0.5, alpha = 0.05, n1 = 50L, n2 = 50L)
  ours <- protocol_text(res, time = as.POSIXct("2026-10-05 12:46:17", tz = "UTC"), provenance = FALSE, eol = "\r\n")
  gpower <- crlf(
    "[1] -- Monday, October 05, 2026 -- 12:46:17",
    "t tests - Means: Difference between two independent means (two groups)",
    "Analysis:\tPost hoc: Compute achieved power ",
    "Input:\tTail(s)\t=\tTwo",
    "\tEffect size d\t=\t0.5",
    "\tα err prob\t=\t0.05",
    "\tSample size group 1\t=\t50",
    "\tSample size group 2\t=\t50",
    "Output:\tNoncentrality parameter δ\t=\t2.5000000",
    "\tCritical t\t=\t1.9844675",
    "\tDf\t=\t98",
    "\tPower (1-β err prob)\t=\t0.6968934")
  expect_identical(strip_end(ours), strip_end(gpower))
})

test_that("two calculations give two numbered records as in G*Power", {
  # GPW-UT-PROT-UTL-U001-UO01--r02, raw protocol SHA-256 9943b6d9...d452f072e.
  r1 <- t2_res("post_hoc", d = 0.5, alpha = 0.05, n1 = 50L, n2 = 50L)
  r2 <- t2_res("post_hoc", d = 0.5, alpha = 0.05, n1 = 60L, n2 = 60L)
  ours <- protocol_text(r1, r2, time = as.POSIXct(c("2026-10-05 13:57:49", "2026-10-05 13:57:50"), tz = "UTC"),
                        provenance = FALSE, eol = "\r\n")
  lines <- strsplit(ours, "\r\n", fixed = TRUE)[[1]]
  expect_identical(lines[1], "[1] -- Monday, October 05, 2026 -- 13:57:49")
  expect_identical(lines[13], "")
  expect_identical(lines[14], "[2] -- Monday, October 05, 2026 -- 13:57:50")
  expect_identical(lines[20:26],
                   c("\tSample size group 1\t=\t60", "\tSample size group 2\t=\t60",
                     "Output:\tNoncentrality parameter δ\t=\t2.7386128", "\tCritical t\t=\t1.9802722",
                     "\tDf\t=\t118", "\tPower (1-β err prob)\t=\t0.7752659", ""))
})

test_that("an a priori record matches the G*Power protocol export", {
  # GPW-T-MEAN2-APR-B001-O01--r01, raw protocol SHA-256 e229f324...d9cd4ae.
  res <- t2_res("a_priori", tails = "one", d = 0.5, alpha = 0.05, power = 0.95, allocation_ratio = 1)
  ours <- protocol_text(res, time = as.POSIXct("2026-10-05 12:18:11", tz = "UTC"), provenance = FALSE, eol = "\r\n")
  gpower <- crlf(
    "[1] -- Monday, October 05, 2026 -- 12:18:11",
    "t tests - Means: Difference between two independent means (two groups)",
    "Analysis:\tA priori: Compute required sample size ",
    "Input:\tTail(s)\t=\tOne",
    "\tEffect size d\t=\t0.5",
    "\tα err prob\t=\t0.05",
    "\tPower (1-β err prob)\t=\t0.95",
    "\tAllocation ratio N2/N1\t=\t1",
    "Output:\tNoncentrality parameter δ\t=\t3.3166248",
    "\tCritical t\t=\t1.6536580",
    "\tDf\t=\t174",
    "\tSample size group 1\t=\t88",
    "\tSample size group 2\t=\t88",
    "\tTotal sample size\t=\t176",
    "\tActual power\t=\t0.9514254")
  expect_identical(strip_end(ours), strip_end(gpower))
})

test_that("records of the other analysis types carry their live labels", {
  rec <- protocol_record(t2_res("sensitivity", alpha = 0.05, power = 0.8, n1 = 50L, n2 = 50L))
  expect_identical(rec$analysis_line, "Sensitivity: Compute required effect size ")
  expect_identical(rec$outputs$label, c("Noncentrality parameter δ", "Critical t", "Df", "Effect size d"))
  rec <- protocol_record(t2_res("criterion", d = 0.5, power = 0.8, n1 = 50L, n2 = 50L))
  expect_identical(rec$analysis_line, "Criterion: Compute required α ")
  expect_identical(rec$outputs$label[4], "α err prob")
  rec <- protocol_record(t2_res("compromise", d = 0.5, error_ratio = 1, n1 = 50L, n2 = 50L))
  expect_identical(rec$analysis_line, "Compromise: Compute implied α & power ")
  expect_identical(rec$inputs$label, c("Tail(s)", "Effect size d", "β/α ratio", "Sample size group 1", "Sample size group 2"))
  expect_identical(rec$outputs$label[4:6], c("α err prob", "β err prob", "Power (1-β err prob)"))
  expect_length(rec$unresolved, 0L)
})

test_that("numbers are formatted as the G*Power output fields display them", {
  # Rule derived from the 39,523 real-valued displays of the emulator corpus
  # (method notes); examples are displays of the cited records.
  f <- inferguidR:::.ig_gpower_number
  expect_identical(f(0.98290120), "0.9829012")                 # GPW-EX-UNCOFF-CRI-B124-D01--r01
  expect_identical(f(-1.959964), "-1.9599640")                  # GPW-Z-POIS-POS-B121-O02--r01
  expect_identical(f(954511.7), "954511.7")                     # Denominator df display
  expect_identical(f(1665.318), "1665.318")                     # GPW-T-WMW-APR-B107-O01--r01
  expect_identical(f(0.001046192), "0.001046192")               # GPW-EX-UNCOND-APR-B046-D01--r01
  expect_identical(f(9.983523e-05), "9.983523e-05")             # Actual alpha display
  expect_identical(f(2.511521e9), "2.511521e+09")               # GPW-T-WSRP-COM-B106-D01--r01
  expect_identical(f(0.0466560), "0.0466560")                   # GPW-EX-UNCOFF-CRI-B259-E01--r01
  expect_identical(f(0), "0")
  expect_identical(f(Inf), "inf")
  expect_identical(f(NaN), "-nan(ind)")
})

test_that("family and analysis texts cover the registry", {
  fam <- inferguidR:::.ig_protocol_family_text
  meta_fam <- unique(vapply(inferguidR:::.ig_list_method_ids(), function(id) inferguidR:::.ig_get_metadata(id)$family, ""))
  expect_true(all(setdiff(meta_fam, "boschloo") %in% names(fam)))
  expect_identical(inferguidR:::.ig_protocol_analysis_text("t_generic", "post_hoc"), "Post hoc: Compute power ")
  expect_identical(inferguidR:::.ig_protocol_analysis_text("f_generic", "sensitivity"), "Sensitivity: Compute noncentrality parameter ")
  expect_identical(inferguidR:::.ig_protocol_analysis_text(t2, "compromise"), "Compromise: Compute implied α & power ")
  # specs/gpower-app-capture.json is a development file; when present, every
  # procedure's live texts must agree with the tables.
  cap <- file.path(th_pkg_root(), "specs", "gpower-app-capture.json")
  skip_if_not(file.exists(cap), "development specs not available")
  procs <- jsonlite::read_json(cap)$procedures
  for (id in names(procs)) {
    p <- procs[[id]]
    meta <- inferguidR:::.ig_get_metadata(id)
    expect_identical(unname(fam[meta$family]), p$test_family_live, info = id)
    expect_identical(meta$label_raw$windows_3197_live, p$statistical_test_label_live, info = id)
    for (a in unlist(p$analysis_types_live)) {
      mode <- c("A priori" = "a_priori", "Post hoc" = "post_hoc", "Sensitivity" = "sensitivity",
                "Criterion" = "criterion", "Compromise" = "compromise")[[sub(":.*$", "", a)]]
      expect_identical(inferguidR:::.ig_protocol_analysis_type(id, mode), a, info = paste(id, mode))
    }
  }
})

test_that("labels without a value are reported, never invented", {
  fake <- list(quantities = list(ncp = 2, critical = c(-1.9, 1.9), df = 10), achieved = list(power = 0.6),
               solved = list(power = 0.6), group_sizes = NULL, method_outputs = list())
  out <- inferguidR:::.ig_protocol_generic_outputs(fake, c("Noncentrality parameter λ", "Critical F", "Df",
                                                           "Power (1-β err prob)", "Some other output"), NULL)
  expect_identical(names(out), c("Noncentrality parameter λ", "Critical F", "Df", "Power (1-β err prob)"))
  expect_identical(out[["Critical F"]], 1.9)
})

test_that("a result without a number gets a record without outputs", {
  # d = 1e-9 cannot reach power 0.95 within the admissible designs: infeasible.
  res <- t2_res("a_priori", d = 1e-9, alpha = 0.05, power = 0.95, allocation_ratio = 1)
  skip_if(identical(res$status$code, "ok"), "the engine solves this case")
  txt <- protocol_text(res, time = as.POSIXct("2026-01-01 00:00:00", tz = "UTC"), provenance = FALSE)
  expect_match(txt, paste0("Output:\t\\(no result: ", res$status$code))
  expect_false(grepl("Total sample size", txt, fixed = TRUE))
})

test_that("provenance, options and argument checks", {
  res <- t2_res("post_hoc", d = 0.5, alpha = 0.05, n1 = 50L, n2 = 50L)
  txt <- protocol_text(res, time = as.POSIXct("2026-10-05 12:00:00", tz = "UTC"))
  expect_match(txt, "inferguidR:\tversion ")
  expect_match(txt, "method status independently_checked|method status gpower_validated|method status development")
  expect_match(txt, res$spec_hash, fixed = TRUE)
  expect_match(txt, "inferguidR::calculate_power(", fixed = TRUE)
  expect_false(grepl("\r", txt, fixed = TRUE))
  txt2 <- protocol_text(res, res, time = as.POSIXct("2026-10-05 12:00:00", tz = "UTC"), start = 7L, provenance = FALSE)
  expect_match(txt2, "^\\[7\\] -- ")
  expect_match(txt2, "\n\\[8\\] -- ")
  rec <- protocol_record(res)
  expect_s3_class(rec, "ig_protocol_record")
  expect_null(rec$options_line)
  expect_identical(protocol_text(rec, time = as.POSIXct("2026-10-05 12:00:00", tz = "UTC"), provenance = FALSE),
                   protocol_text(res, time = as.POSIXct("2026-10-05 12:00:00", tz = "UTC"), provenance = FALSE))
  expect_error(protocol_text(), class = "inferguidR_invalid_input")
  expect_error(protocol_text(1), class = "inferguidR_invalid_input")
  expect_error(protocol_text(res, eol = "\r"), class = "inferguidR_invalid_input")
  expect_error(protocol_text(res, time = "now"), class = "inferguidR_invalid_input")
  expect_error(protocol_text(res, start = 0), class = "inferguidR_invalid_input")
  expect_error(protocol_record(list()), class = "inferguidR_invalid_input")
})

test_that("dates are written with English names independent of the locale", {
  # The header uses fixed English weekday and month names (not strftime), so
  # the locale cannot change it; the C locale gives the reference text.
  old <- Sys.getlocale("LC_TIME")
  on.exit(Sys.setlocale("LC_TIME", old), add = TRUE)
  Sys.setlocale("LC_TIME", "C")
  for (tm in as.POSIXct(c("2026-03-01 09:05:07", "2026-10-05 13:57:49", "2027-12-31 23:59:59"), tz = "UTC")) {
    tm <- as.POSIXct(tm, origin = "1970-01-01", tz = "UTC")
    expect_identical(inferguidR:::.ig_protocol_header(3, tm),
                     paste0("[3] -- ", format(tm, "%A, %B %d, %Y -- %H:%M:%S", tz = "UTC")))
  }
})

test_that("a live label that shortens a parameter label is filled from that parameter", {
  id <- "t_means_matched_pairs"
  skip_if(is.null(inferguidR:::.ig_get_method(id)), "no engine for the matched-pairs t test")
  res <- tryCatch(calculate_power(analysis_spec(id, "post_hoc", inputs = list(dz = 0.5, alpha = 0.05, n_total = 30L),
                                                options = list(tails = "two"))), error = function(e) NULL)
  skip_if(is.null(res) || !identical(res$status$code, "ok"), "the matched-pairs specification differs")
  rec <- protocol_record(res)
  # Live label "Total sample size"; parameter label "Total sample size (number of pairs)".
  expect_identical(rec$inputs$label, c("Tail(s)", "Effect size dz", "α err prob", "Total sample size"))
  expect_identical(rec$inputs$text[4], "30")
  expect_length(rec$unresolved, 0L)
})

test_that("a generic 'Options:' line says that it is not G*Power's wording", {
  meta <- list(options = list(
    list(option_id = "tails", kind = "main_window_input", value_live = list(two = list(text = "Two"))),
    list(option_id = "method", kind = "options_dialog", ui_label = "Method",
         value_live = list(a = list(text = "Wilks U                           (Rao, 1951)"))),
    list(option_id = "limit", kind = "options_dialog", ui_label = "Limit", value_live = list())))
  spec <- list(options = list(tails = "two", method = "a", limit = 1000))
  expect_identical(inferguidR:::.ig_protocol_generic_options(spec, meta), "Wilks U (Rao, 1951), Limit = 1000")
  expect_null(inferguidR:::.ig_protocol_generic_options(list(options = list(tails = "two")), meta))
  # A method without dialog options has no 'Options:' line and no such note.
  rec <- protocol_record(t2_res("post_hoc", d = 0.5, alpha = 0.05, n1 = 50L, n2 = 50L))
  expect_null(rec$options_line)
  expect_false(any(grepl("Options:", rec$notes, fixed = TRUE)))
  # A method with dialog options gets the line and the note (G*Power's own
  # text for this state is 'α balancing: α/2 on each side', GPW-EX-BINOM records).
  id <- "exact_proportion_one_sample_binomial"
  skip_if(is.null(inferguidR:::.ig_get_method(id)), "no engine for the binomial test")
  res <- tryCatch(calculate_power(analysis_spec(id, "post_hoc", inputs = list(g = 0.15, p0 = 0.5, alpha = 0.05, n_total = 40L),
                                                options = list(tails = "two", alpha_balancing = "half_each_side"))),
                  error = function(e) NULL)
  skip_if(is.null(res) || !identical(res$status$code, "ok"), "the binomial specification differs")
  rec <- protocol_record(res)
  expect_false(is.null(rec$options_line))
  expect_true(any(grepl("G*Power writes its own wording", rec$notes, fixed = TRUE)))
  txt <- protocol_text(rec, time = as.POSIXct("2026-10-05 12:00:00", tz = "UTC"))
  expect_match(txt, "\nOptions:\t")
  expect_match(txt, "note: the 'Options:' line lists the chosen options by their dialog texts", fixed = TRUE)
})
