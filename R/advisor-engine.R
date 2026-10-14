# Rule-graph interpreter -------------------------------------------------------
#
# Semantics (rules.json, section `engine`):
#
# * The walk starts at engine.start_question.
# * The effective answer of an option is the option itself, except for a
#   'dont_know' option whose policy is route_without_assumption or
#   route_with_assumption: its effective answer is `routes_as`. The raw answer
#   is always kept.
# * After an answer, the transitions whose `from` is the current question are
#   evaluated in ascending `order` against the effective answers given so far;
#   the first transition whose condition holds names the next question or the
#   outcome. Every question ends with an unconditional guard transition to
#   out.needs_information; a valid rule set never reaches a guard.
# * On arrival at a question, its `applies_when` condition must hold (only
#   applicable follow-ups are asked); a violation is a rule-set defect.
# * A walk that stops at an unanswered question yields out.needs_information
#   with that question as the pending question (incomplete_policy).

# Effective answer of one option record: its `routes_as` for a routing
# 'dont_know', otherwise its own option id.
.ig_advisor_effective <- function(option) {
  pol <- option$dont_know_policy
  if (identical(option$option_id, "dont_know") && .ig_is_string(pol) &&
      pol %in% c("route_without_assumption", "route_with_assumption")) {
    return(option$routes_as)
  }
  option$option_id
}

# Index of the first transition from `qid` whose condition holds.
.ig_advisor_first_transition <- function(rs, qid, eff) {
  fs <- rs$compiled$transitions[[qid]]
  for (i in seq_along(fs)) if (fs[[i]](eff)) return(i)
  NA_integer_
}

# Replay an ordered list of raw answers (list of list(question_id, option_id)).
# Returns the walk: steps, effective and raw answer vectors, the node reached
# (question or outcome) and, if the walk stopped at a question, that question
# as `pending`. Answers that do not follow the graph are an invalid_input error.
.ig_advisor_walk <- function(rs, answers) {
  node <- rs$docs$rules$engine$start_question
  eff <- stats::setNames(character(), character())
  raw <- eff
  steps <- list()
  visited <- character()
  i <- 0L
  repeat {
    if (.ig_advisor_is_outcome(rs, node)) break
    q <- rs$questions[[node]]
    if (is.null(q)) {
      .ig_abort("invalid_input", sprintf("rule set defect: transition target '%s' does not exist", node))
    }
    if (node %in% visited) {
      .ig_abort("invalid_input", sprintf("rule set defect: question '%s' is reached twice (cycle)", node))
    }
    if (!rs$compiled$applies[[node]](eff)) {
      .ig_abort("invalid_input", sprintf(
        "rule set defect: question '%s' was reached although its applies_when condition does not hold", node))
    }
    visited <- c(visited, node)
    if (i >= length(answers)) break
    i <- i + 1L
    a <- answers[[i]]
    if (!identical(a$question_id, node)) {
      .ig_abort("invalid_input", sprintf(
        "answer %d is for question '%s', but the walkthrough is at question '%s'",
        i, format(a$question_id), node),
        list(problems = list(sprintf("expected an answer to '%s'", node)),
             expected_question_id = node))
    }
    opt <- rs$options[[node]][[a$option_id]]
    if (is.null(opt)) {
      .ig_abort("invalid_input", sprintf("'%s' is not an option of question '%s'",
                                         format(a$option_id), node))
    }
    e <- .ig_advisor_effective(opt)
    eff[[node]] <- e
    raw[[node]] <- a$option_id
    k <- .ig_advisor_first_transition(rs, node, eff)
    if (is.na(k)) {
      .ig_abort("invalid_input", sprintf("rule set defect: no transition from '%s' applies", node))
    }
    tr <- rs$transitions[[node]][[k]]
    steps[[length(steps) + 1L]] <- list(
      step = i, question_id = node, raw_answer = a$option_id, effective_answer = e,
      dont_know_policy = if (identical(a$option_id, "dont_know")) opt$dont_know_policy else NULL,
      transition_id = tr$transition_id, transition_order = as.integer(tr$order),
      kind = tr$kind, to = tr$to
    )
    node <- tr$to
  }
  if (i < length(answers)) {
    extra <- vapply(answers[(i + 1L):length(answers)], function(a) format(a$question_id), "")
    .ig_abort("invalid_input", sprintf(
      "%d answer(s) follow the end of the walkthrough (%s): %s",
      length(answers) - i, node, paste(extra, collapse = ", ")))
  }
  outcome <- .ig_advisor_is_outcome(rs, node)
  list(node = node,
       outcome_id = if (outcome) node else NULL,
       pending = if (outcome) NULL else node,
       steps = steps, effective = eff, raw = raw, visited = visited,
       guard_reached = any(vapply(steps, function(s) identical(s$kind, "guard"), logical(1))))
}

# Status rules (rules.json engine.status_rules), applied to an outcome reached
# by `walk`. `override` is NULL or an allowed override variant id.
.ig_advisor_status <- function(rs, outcome, walk, override = NULL) {
  base <- outcome$status
  status <- base
  applied <- "sr.base"
  downgrade <- function(s) if (identical(s, "recommended")) "conditional" else s
  assumptions <- list()
  dk <- Filter(function(s) identical(s$raw_answer, "dont_know") &&
                 identical(s$dont_know_policy, "route_with_assumption"), walk$steps)
  if (length(dk)) {
    applied <- c(applied, "sr.dont_know_downgrade")
    for (s in dk) {
      opt <- rs$options[[s$question_id]][["dont_know"]]
      assumptions[[length(assumptions) + 1L]] <- list(
        assumption_id = paste0("dk.", s$question_id), state = "unknown",
        text_en = opt$assumption_en, text_de = opt$assumption_de,
        source = "dont_know", question_id = s$question_id
      )
      if (!(status %in% c("needs_information", "outside_supported_scope"))) status <- downgrade(status)
    }
  }
  mods <- character()
  conds <- rs$compiled$modifiers[[outcome$outcome_id]]
  for (j in seq_along(outcome$status_modifiers)) {
    m <- outcome$status_modifiers[[j]]
    if (conds[[j]](walk$effective)) {
      mods <- c(mods, m$modifier_id)
      assumptions[[length(assumptions) + 1L]] <- list(
        assumption_id = m$modifier_id, state = "unknown",
        text_en = m$assumption_en, text_de = m$assumption_de,
        source = "status_modifier", question_id = NULL
      )
      if (identical(m$effect, "conditional") &&
          !(status %in% c("needs_information", "outside_supported_scope"))) {
        status <- downgrade(status)
      }
    }
  }
  if (length(mods)) applied <- c(applied, "sr.modifiers")
  if (!is.null(override)) {
    applied <- c(applied, "override_policy")
    status <- downgrade(status)
  }
  applied <- c(applied, "sr.no_upgrade")
  list(base_status = base, status = status, rules_applied = applied,
       modifiers_applied = mods, assumptions_made = assumptions)
}

# Outcome record for a walk (the pending case uses out.needs_information).
.ig_advisor_outcome_of <- function(rs, walk) {
  if (!is.null(walk$outcome_id)) return(rs$outcomes[[walk$outcome_id]])
  rs$outcomes[["out.needs_information"]]
}
