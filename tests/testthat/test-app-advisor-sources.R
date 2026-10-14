test_that("offline source trees expose every audited node without changing routing", {
 skip_if_not_installed("shiny")
 rules <- advisor_rules()
 r <- advisor_recommendation(advisor_session())
 before <- serialize(r,NULL)
 for (lang in c("en","de")) {
  html <- as.character(.ig_advisor_source_view(r,lang,rules))
  nodes <- c(rules$docs$uzh_tree$nodes,rules$docs$descrtab2_tree$nodes,rules$docs$descrtab2_tree$code_only_nodes)
  for (node in nodes) expect_match(html,paste0('data-source-node="',node$node_id %||% node$id,'"'),fixed=TRUE)
  expect_match(html,"https://",fixed=TRUE)
  expect_false(grepl("<iframe",html,fixed=TRUE))
 }
 expect_identical(serialize(r,NULL),before)
})

test_that("source comparison retains original and selected methods after override", {
 skip_if_not_installed("shiny")
 map <- c(q.aim="compare",q.stage="analysis",q.outcome_type="continuous",q.unit="individuals",q.structure="independent",
   q.independence_check="no",q.n_factors="one",q.n_groups="two",q.covariates="no",q.estimand_continuous="mean",
   q.variance="not_justified",q.direction="two_sided",q.missingness="none")
 state <- advisor_replay(map)
 rules <- advisor_rules()
 for (override in c(FALSE,TRUE)) {
  s <- if(override) advisor_override(state,"rank_sum_mann_whitney") else state
  r <- advisor_recommendation(s)
  html <- as.character(.ig_advisor_source_view(r,"en",rules))
  expect_match(html,'data-crosswalk-id="cw002"',fixed=TRUE)
  expect_match(html,"without a Levene gate",fixed=TRUE)
  expect_match(html,"not a replay of your answers",fixed=TRUE)
  if(override) {
   expect_match(html,'data-crosswalk-id="cw004"',fixed=TRUE)
   expect_match(html,"Original rule recommendation",fixed=TRUE)
   expect_match(html,"Selected method",fixed=TRUE)
  } else expect_false(grepl('data-crosswalk-id="cw004"',html,fixed=TRUE))
 }
})

test_that("source views render every acceptance outcome in both languages", {
 skip_if_not_installed("shiny")
 rules <- advisor_rules()
 scenarios <- jsonlite::fromJSON(system.file("advisor","acceptance-scenarios.json",package="inferguidR"),simplifyVector=FALSE)$scenarios
 for(sc in scenarios) {
  if(!length(sc$answers)) next
  answers <- stats::setNames(vapply(sc$answers,function(a) a$answer,""),vapply(sc$answers,function(a) a$question_id,""))
  r <- advisor_recommendation(advisor_replay(answers))
  for(lang in c("en","de")) {
   expect_no_error(html <- as.character(.ig_advisor_source_view(r,lang,rules)))
   expect_match(html,"ig-source-comparison",fixed=TRUE)
  }
 }
})

test_that("incomplete outcomes do not select unrelated empty-method crosswalk rows", {
 skip_if_not_installed("shiny")
 html <- as.character(.ig_advisor_source_view(advisor_recommendation(advisor_session()),"en",advisor_rules()))
 expect_match(html,'data-crosswalk-id="cw042"',fixed=TRUE)
 expect_false(grepl('data-crosswalk-id="cw044"',html,fixed=TRUE))
 expect_false(grepl('data-crosswalk-id="cw045"',html,fixed=TRUE))
 expect_match(html,"No method is specified by this entry",fixed=TRUE)
 expect_false(grepl("<h4></h4>",html,fixed=TRUE))
})
