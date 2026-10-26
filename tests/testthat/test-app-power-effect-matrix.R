test_that("every advertised helper layout preserves locale apply and stale contracts", {
  skip_if_not_installed("shiny")
  samples <- list(mean_1 = 2.5, mean_2 = 1.25, sd = 2.5, sd_1 = 2.5, sd_2 = 2,
    r = .3, rho = .3, n1 = 20, n2 = 24, mean_h0 = 0, mean_h1 = 1.25,
    mean_diff = 1.25, sd_diff = 2.5, var0 = 1.25, var1 = 2.5,
    means = c(1.25, 2.5, 3.75), sizes = c(20, 24, 28),
    var_effect = .5, var_error = 2.5, partial_eta_squared = .2)
  covered <- 0L
  for (method in list_methods()$id) {
    layouts <- .ig_power_effect_layouts(method)
    for (input_mode in names(layouts)) {
      entered <- samples[layouts[[input_mode]]]
      expected <- effect_size(method, input_mode, inputs = entered)
      form_value <- shiny::reactiveVal(.ig_workbench_form(method, "post_hoc"))
      locale <- shiny::reactiveVal("en")
      applied <- new.env(parent = emptyenv())
      applied$value <- NULL
      input_names <- paste(method, input_mode, names(entered), sep = "_")
      english <- stats::setNames(lapply(entered, function(v)
        paste(vapply(v, .ig_format_double, ""), collapse = ";")), input_names)
      german <- lapply(english, function(v) gsub(".", ",", v, fixed = TRUE))
      shiny::testServer(.ig_power_effect_server,
        args = list(form = form_value, language = shiny::reactive(locale()),
          apply_value = function(field, value) applied$value <- list(field = field, value = value)), {
        session$flushReact()
        session$setInputs(mode = input_mode)
        expect_no_error(output$controls)
        expect_no_error(output$fields)
        do.call(session$setInputs, c(english, list(calculate = 1)))
        expect_identical(current()$status$code, "ok", info = paste(method, input_mode))
        expect_equal(current()$value, expected$value)
        expect_no_error(output$result)
        apply_id <- paste0('id="', session$ns("apply"), '"')
        expect_match(output$actions$html, apply_id, fixed = TRUE)
        expect_null(applied$value)
        session$setInputs(apply = 1)
        expect_identical(applied$value, list(field = expected$effect, value = expected$value))
        applied$value <- NULL
        do.call(session$setInputs, c(stats::setNames(list(""), input_names[1]), list(apply = 2)))
        expect_null(current())
        expect_null(applied$value)
        expect_match(output$message, "current, applicable", fixed = TRUE)
        expect_false(grepl(apply_id, output$actions$html, fixed = TRUE))
        do.call(session$setInputs, c(english, list(calculate = 2)))
        locale("de")
        session$flushReact()
        rendered <- output$fields$html
        for (value in german) expect_match(rendered, paste0('value="', value, '"'), fixed = TRUE)
        # A browser sends the localized values after renderUI replaces the controls.
        do.call(session$setInputs, german)
        expect_equal(current()$value, expected$value)
        expect_match(output$actions$html, "In Entwurf uebernehmen", fixed = TRUE)
        expect_null(applied$value)
        session$setInputs(apply = 3)
        expect_identical(applied$value, list(field = expected$effect, value = expected$value))
      })
      covered <- covered + 1L
    }
  }
  expect_equal(covered, 20L)
})
