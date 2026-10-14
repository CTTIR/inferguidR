test_that("effect helper parses locale numbers without evaluating code", {
 expect_equal(.ig_power_effect_numbers("0,2;0,5", "de", TRUE),c(.2,.5))
 expect_equal(.ig_power_effect_numbers("-2.5", "en"),-2.5)
 for(x in c("", "0.5;", "0.5;0.8", "system('true')", "Inf"))
  expect_error(.ig_power_effect_numbers(x,"en"),class="inferguidR_invalid_input")
})

test_that("effect helpers preserve signed values and never commit a result", {
 skip_if_not_installed("shiny")
 applied <- NULL
 f <- shiny::reactiveVal(.ig_workbench_form("t_means_one_sample","post_hoc",list(tails="two")))
 shiny::testServer(.ig_power_effect_server,args=list(form=f,language=shiny::reactive("en"),
   apply_value=function(field,value) applied <<- list(field=field,value=value)), {
  session$flushReact()
  session$setInputs(mode="from_means",t_means_one_sample_from_means_mean_h0="10",
    t_means_one_sample_from_means_mean_h1="8",t_means_one_sample_from_means_sd="4",calculate=1)
  expect_equal(current()$value,-.5)
  expect_null(applied)
  session$setInputs(apply=1)
  expect_identical(applied,list(field="d",value=-.5))
  applied <- NULL
  session$setInputs(t_means_one_sample_from_means_sd="0",apply=2)
  expect_null(current())
  expect_null(applied)
  session$setInputs(calculate=2)
  expect_match(message(),"positive")
  session$setInputs(t_means_one_sample_from_means_sd="4",calculate=3)
  f(.ig_workbench_form("t_means_one_sample","sensitivity",list(tails="two")))
  session$flushReact()
  session$setInputs(apply=3)
  expect_null(applied)
  expect_match(output$message,"applicable")
 })
})

test_that("all helper layouts call their exact supported effect API", {
 samples <- list(mean_1=2,mean_2=1,sd=2,sd_1=2,sd_2=2,r=.3,rho=.3,n1=20,n2=20,
  mean_h0=0,mean_h1=1,mean_diff=1,sd_diff=2,var0=1,var1=2,means=c(1,2,3),sizes=c(20,20,20),
  var_effect=.5,var_error=2,partial_eta_squared=.2)
 ids <- list_methods()$id
 covered <- 0L
 for(id in ids) {
  layouts <- .ig_power_effect_layouts(id)
  for(mode in names(layouts)) {
   r <- effect_size(id,mode,inputs=samples[layouts[[mode]]])
   expect_identical(r$status$code,"ok",info=paste(id,mode))
   expect_true(is.finite(r$value),info=paste(id,mode))
   covered <- covered+1L
  }
 }
 expect_equal(covered,20)
})

test_that("effect helper exposes numerical failure and never applies it", {
 skip_if_not_installed("shiny")
 applied <- NULL
 f <- shiny::reactiveVal(.ig_workbench_form("f_variances_two_groups","post_hoc"))
 shiny::testServer(.ig_power_effect_server,args=list(form=f,language=shiny::reactive("en"),
   apply_value=function(field,value) applied <<- list(field=field,value=value)), {
  session$flushReact()
  session$setInputs(mode="from_variances",f_variances_two_groups_from_variances_var0="1e-308",
    f_variances_two_groups_from_variances_var1="1e308",calculate=1)
  expect_identical(current()$status$code,"numerical_failure")
  expect_match(output$result$html,"numerical_failure")
  expect_match(output$result$html,"no finite value")
  session$setInputs(apply=1)
  expect_null(applied)
 })
})
