test_that("dynamic power inputs preserve drafts and their entry locale while unbound", {
 skip_if_not_installed("shiny");skip_if_not_installed("ggplot2")
 lang<-shiny::reactiveVal("en")
 shiny::testServer(.ig_power_workspace_server,args=list(language=lang),{
  session$flushReact()
  session$setInputs(method="boschloo_proportions_two_independent",mode="post_hoc",
   option_boschloo_proportions_two_independent_variant="central",research=TRUE,
   field_boschloo_proportions_two_independent_p1="0.2",
   field_boschloo_proportions_two_independent_p2="0.7",
   field_boschloo_proportions_two_independent_alpha="0.05",
   field_boschloo_proportions_two_independent_n1="3",
   field_boschloo_proportions_two_independent_n2="5")
  session$setInputs(calculate=1)
  expect_identical(current()$status$code,"ok")
  layout<-field_layout()
  option_editor<-option_layout();target_editor<-target_layout()
  session$setInputs(option_boschloo_proportions_two_independent_variant="minlike")
  expect_identical(field_layout(),layout)
  expect_identical(option_layout(),option_editor)
  expect_identical(target_layout(),target_editor)
  session$setInputs(option_boschloo_proportions_two_independent_variant="central")
  expect_identical(field_layout(),layout)
  expect_identical(option_layout(),option_editor)
  original<-state()$committed$analysis_hash
  session$setInputs(field_boschloo_proportions_two_independent_p1=NULL,
    field_boschloo_proportions_two_independent_n1=NULL)
  lang("de");session$flushReact()
  expect_match(output$fields$html,'value="0,2"',fixed=TRUE)
  expect_match(output$fields$html,'value="3"',fixed=TRUE)
  expect_identical(state()$committed$analysis_hash,original)
  # Until the browser rebinds matching values the old result is still stale.
  expect_true(state()$dirty)
  session$setInputs(field_boschloo_proportions_two_independent_p1="invalid",
    field_boschloo_proportions_two_independent_n1="")
  session$setInputs(field_boschloo_proportions_two_independent_p1=NULL,
    field_boschloo_proportions_two_independent_n1=NULL)
  lang("en");session$flushReact()
  expect_match(output$fields$html,'value="invalid"',fixed=TRUE)
  expect_identical(field_drafts$field_boschloo_proportions_two_independent_n1$text,"")
  expect_true(state()$dirty)
  session$setInputs(method="t_means_one_sample")
  session$setInputs(method="boschloo_proportions_two_independent")
  expect_match(output$fields$html,'value="invalid"',fixed=TRUE)
  expect_identical(state()$committed$analysis_hash,original)
 })
})

test_that("effect helper tolerates unsupported methods and returns to supported ones", {
 skip_if_not_installed("shiny")
 f<-shiny::reactiveVal(.ig_workbench_form("t_means_one_sample","post_hoc",list(tails="two")))
 shiny::testServer(.ig_power_effect_server,args=list(form=f,language=shiny::reactive("en"),apply_value=function(...)NULL),{
  session$flushReact();expect_length(fields(),3)
  expect_type(output$fields,"list")
  for(id in c("boschloo_proportions_two_independent","f_generic")) {
   f(.ig_workbench_form(id,"post_hoc"));session$flushReact()
   expect_identical(fields(),character())
   expect_no_error(output$fields);expect_no_error(output$actions);expect_no_error(output$result)
   expect_match(output$controls$html,"does not yet have an input helper",fixed=TRUE)
  }
  f(.ig_workbench_form("t_means_one_sample","post_hoc",list(tails="two")));session$flushReact()
  expect_length(fields(),3);expect_type(output$fields,"list")
 })
})
