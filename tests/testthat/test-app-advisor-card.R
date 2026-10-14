test_that("recommendation cards render acceptance scenarios in both languages", {
 skip_if_not_installed("shiny")
 scenarios <- jsonlite::fromJSON(system.file("advisor","acceptance-scenarios.json",package="inferguidR"),simplifyVector=FALSE)$scenarios
 rules <- advisor_rules()
 for (sc in scenarios) {
  if (!length(sc$answers)) next
  answers <- stats::setNames(vapply(sc$answers,function(a) a$answer,""),vapply(sc$answers,function(a) a$question_id,""))
  s <- advisor_replay(answers)
  r <- advisor_recommendation(s)
  for (lang in c("en","de")) {
   html <- as.character(.ig_advisor_card(r,lang,rules))
   expect_match(html,"ig-recommendation",info=paste(sc$scenario_id,lang))
   expect_match(html,r$rules$rules_version,fixed=TRUE)
   if (!is.null(r$preferred)) expect_true(grepl(r$preferred[[paste0("name_",lang)]],html,fixed=TRUE),info=sc$scenario_id)
  }
 }
})

test_that("cards preserve resource limitations and safe links", {
 skip_if_not_installed("shiny")
 map <- c(q.aim="compare",q.stage="analysis",q.outcome_type="continuous",q.unit="individuals",q.structure="independent",
   q.independence_check="no",q.n_factors="one",q.n_groups="two",q.covariates="no",q.estimand_continuous="mean",
   q.variance="not_justified",q.direction="two_sided",q.missingness="none")
 r <- advisor_recommendation(advisor_replay(map))
 before <- r$recommendation_hash
 html <- as.character(.ig_advisor_card(r,"en",advisor_rules()))
 expect_match(html,"not assessed")
 expect_match(html,"effect size cohens d")
 expect_false(grepl("exact registry match",html,fixed=TRUE))
 expect_match(html,"supporting concepts")
 expect_match(html,"R template is bundled with the package")
 expect_match(html,"never numerical validation evidence")
 expect_match(html,"External methodological background")
 expect_identical(r$recommendation_hash,before)
 expect_false(grepl("href",as.character(.ig_advisor_link("unsafe","javascript:alert(1)"))))
 expect_match(as.character(.ig_advisor_link("<script>","https://example.org")),"&lt;script&gt;",fixed=TRUE)
})

test_that("the guided app applies only rule-allowed alternatives", {
 skip_if_not_installed("shiny")
 map <- c(q.aim="compare",q.stage="analysis",q.outcome_type="continuous",q.unit="individuals",q.structure="independent",
   q.independence_check="no",q.n_factors="one",q.n_groups="two",q.covariates="no",q.estimand_continuous="mean",
   q.variance="not_justified",q.direction="two_sided",q.missingness="none")
 shiny::testServer(.ig_advisor_server,args=list(language=shiny::reactive("en")), {
  state(advisor_replay(map));session$flushReact()
  expect_match(output$recommendation$html,"Assumptions to assess")
  expect_match(output$override_controls$html,"rank_sum_mann_whitney")
  session$setInputs(override_variant="rank_sum_mann_whitney",apply_override=1)
  r <- advisor_recommendation(state())
  expect_identical(r$status,"conditional")
  expect_identical(r$preferred$analysis_variant_id,"rank_sum_mann_whitney")
  expect_match(output$recommendation$html,"Conditional recommendation")
  before <- state()
  session$setInputs(override_variant="made_up",apply_override=2)
  expect_identical(state(),before)
  expect_true(nzchar(message()))
 })
})

test_that("override cards keep original and selected assumption provenance distinct", {
 skip_if_not_installed("shiny")
 map <- c(q.aim="compare",q.stage="analysis",q.outcome_type="continuous",q.unit="individuals",q.structure="independent",
   q.independence_check="no",q.n_factors="one",q.n_groups="two",q.covariates="no",q.estimand_continuous="mean",
   q.variance="not_justified",q.direction="two_sided",q.missingness="none")
 state <- advisor_replay(map)
 state <- advisor_override(state, "rank_sum_mann_whitney")
 r <- advisor_recommendation(state)
 before <- serialize(r, NULL)
 rules <- advisor_rules()
 for (lang in c("en", "de")) {
  html <- as.character(.ig_advisor_card(r,lang,rules))
  items <- regmatches(html, gregexpr("(?s)<li>.*?</li>", html, perl=TRUE))[[1]]
  for (a in r$assumptions) {
   value <- a[[paste0("text_",lang)]] %||% a$text_en
   item <- items[vapply(items, function(x) grepl(gsub("^<span>|</span>$", "", as.character(shiny::tags$span(value))),x,fixed=TRUE), logical(1))]
   expect_true(length(item) >= 1L)
   rendered <- paste(item,collapse=" ")
   if (identical(a$stated_for,r$preferred$outcome_variant_id))
    expect_match(rendered,if(lang=="en") "Original recommendation:" else "Urspruengliche Empfehlung:",fixed=TRUE)
   if (identical(a$stated_for,r$preferred$analysis_variant_id))
    expect_match(rendered,if(lang=="en") "Selected alternative:" else "Gewaehlte Alternative:",fixed=TRUE)
  }
  expect_match(html,if(lang=="en") "Sources of the original rule recommendation" else "Quellen der urspruenglichen Regelempfehlung",fixed=TRUE)
  expect_match(html,if(lang=="en") "They do not establish support for the selected alternative." else "Sie belegen nicht die gewaehlte Alternative.",fixed=TRUE)
 }
 expect_identical(serialize(r,NULL),before)
 original <- as.character(.ig_advisor_card(advisor_recommendation(advisor_replay(map)),"en",rules))
 expect_match(original,"Methodological sources",fixed=TRUE)
 expect_false(grepl("Sources of the original rule recommendation",original,fixed=TRUE))
})
