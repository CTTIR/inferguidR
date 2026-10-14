test_that("bundled lessons have safe offline paths and matching content hashes", {
  for (lang in c("en", "de")) {
    x <- learning_lesson("boschloo", lang)
    expect_identical(x$analysis_variant_id, "boschloo_two_sided_minlike")
    expect_identical(x$locale, lang)
    for (name in c("markdown", "html", "script")) {
      expect_true(file.exists(x[[name]]))
      expect_identical(digest::digest(file = x[[name]], algo = "sha256"), x$sha256[[name]])
    }
    html <- paste(readLines(x$html, encoding = "UTF-8"), collapse = "\n")
    expect_match(html, "0[.,]004180908")
    expect_false(grepl("<script|<iframe|<img|<link", html, ignore.case = TRUE))
    code <- new.env(parent = baseenv())
    expect_silent(sys.source(x$script, code))
    expect_type(code$run_boschloo_lesson, "closure")
  }
  for (bad in list("../boschloo", "", NA_character_, c("boschloo", "boschloo"), 1, NULL))
    expect_error(learning_lesson(bad), "Unknown bundled lesson")
  for (bad in list("fr", "../en", NA_character_, c("en", "de"), 1, NULL))
    expect_error(learning_lesson(locale = bad), "locale")
})

test_that("advisor lesson tracks the exact selected variant and locale", {
  base <- c(q.aim="compare",q.stage="analysis",q.outcome_type="binary",q.unit="individuals",
    q.structure="independent",q.independence_check="no",q.n_factors="one",q.n_groups="two",
    q.covariates="no",q.binary_margins="group_sizes_fixed",q.direction="two_sided",q.missingness="none")
  central <- advisor_replay(base)
  minlike <- advisor_override(central, "boschloo_two_sided_minlike")
  expect_null(.ig_advisor_lesson(advisor_recommendation(central), "en"))
  for (lang in c("en", "de"))
    expect_identical(.ig_advisor_lesson(advisor_recommendation(minlike), lang)$locale, lang)
  base["q.direction"] <- "one_sided_prespecified"
  expect_null(.ig_advisor_lesson(advisor_recommendation(advisor_replay(base)), "en"))
  skip_if_not_installed("shiny")
  lang <- shiny::reactiveVal("en")
  shiny::testServer(.ig_advisor_server, args = list(language = lang), {
    session$flushReact()
    expect_null(selected_lesson())
    state(minlike); session$flushReact()
    expect_identical(selected_lesson()$locale, "en")
    expect_match(output$worked_lesson$html, "Original worked lesson", fixed = TRUE)
    lang("de"); session$flushReact()
    expect_identical(selected_lesson()$locale, "de")
    expect_match(output$worked_lesson$html, "Eigenstaendige Lektion", fixed = TRUE)
    state(central); session$flushReact()
    expect_null(selected_lesson())
    expect_null(output$worked_lesson)
  })
})
