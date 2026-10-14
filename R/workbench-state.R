# Pure state transitions for the linked workbench. No renderer owns statistics.
# Tickets bind a request to its session, sequence and complete specification.
.ig_workbench_ticket <- function(state) {
  list(session_id = state$session_id, sequence = state$sequence,
       spec_hash = state$committed$spec_hash)
}

.ig_workbench_new <- function(spec, session_id, locale = "en", history_limit = 20L) {
  if (!.ig_is_string(session_id) || !nzchar(session_id))
    .ig_abort("invalid_input", "the workbench needs a nonempty session identifier")
  if (!.ig_is_string(locale) || !locale %in% c("en", "de"))
    .ig_abort("invalid_input", "workbench language must be en or de")
  if (!is.numeric(history_limit) || length(history_limit) != 1L || is.na(history_limit) ||
      !is.finite(history_limit) || history_limit != floor(history_limit) || history_limit < 1 || history_limit > 100)
    .ig_abort("invalid_input", "history_limit must be a whole number from 1 to 100")
  env <- analysis_envelope(spec)
  state <- structure(list(session_id = session_id, sequence = 1, locale = locale,
    committed = env, draft = list(), dirty = FALSE, preview = NULL,
    pending = NULL, history = list(env), cursor = 1L, history_limit = as.integer(history_limit)),
    class = "ig_workbench_state")
  state$pending <- .ig_workbench_ticket(state)
  state
}

.ig_workbench_edit <- function(state, fields) {
  if (!is.list(fields) || !length(fields) || !.ig_is_named(fields) || anyDuplicated(names(fields)))
    .ig_abort("invalid_input", "draft changes must be a named list of fields")
  for (name in names(fields)) state$draft[name] <- fields[name]
  state$dirty <- TRUE
  state
}

.ig_workbench_commit <- function(state, spec) {
  env <- analysis_envelope(spec)
  state$sequence <- state$sequence + 1
  state$committed <- env
  state$draft <- list()
  state$dirty <- FALSE
  state$preview <- NULL
  history <- c(state$history[seq_len(state$cursor)], list(env))
  if (length(history) > state$history_limit) history <- utils::tail(history, state$history_limit)
  state$history <- history
  state$cursor <- length(history)
  state$pending <- .ig_workbench_ticket(state)
  state
}

.ig_workbench_accept <- function(state, ticket, result) {
  if (is.null(state$pending) || !identical(ticket, state$pending) ||
      !identical(result$spec_hash, state$committed$spec_hash)) return(state)
  # The constructor checks the result/spec pairing and existing integrity rules.
  env <- analysis_envelope(state$committed$spec, result)
  state$committed <- env
  state$history[[state$cursor]] <- env
  state$pending <- NULL
  state
}

.ig_workbench_preview <- function(state, envelope) {
  if (!inherits(envelope, "analysis_envelope") || is.null(envelope$result))
    .ig_abort("invalid_input", "a preview needs a complete specification and computed result")
  if (!identical(envelope$analysis_kind, state$committed$analysis_kind))
    .ig_abort("invalid_input", "preview and current workspace must have the same analysis kind")
  state$preview <- analysis_envelope(envelope$spec, envelope$result)
  state
}

.ig_workbench_apply_preview <- function(state) {
  preview <- state$preview
  if (is.null(preview) || !identical(preview$result$status$code, "ok"))
    .ig_abort("invalid_input", "only a successfully computed preview can be applied")
  state <- .ig_workbench_commit(state, preview$spec)
  .ig_workbench_accept(state, state$pending, preview$result)
}

.ig_workbench_history <- function(state, direction) {
  if (!is.numeric(direction) || length(direction) != 1L || is.na(direction) || !direction %in% c(-1, 1))
    .ig_abort("invalid_input", "history direction must be -1 or 1")
  next_cursor <- state$cursor + direction
  if (next_cursor < 1 || next_cursor > length(state$history)) return(state)
  state$sequence <- state$sequence + 1
  state$cursor <- next_cursor
  state$committed <- state$history[[next_cursor]]
  state$draft <- list()
  state$dirty <- FALSE
  state$preview <- NULL
  state$pending <- if (is.null(state$committed$result)) .ig_workbench_ticket(state) else NULL
  state
}

.ig_workbench_locale <- function(state, locale) {
  if (!.ig_is_string(locale) || !locale %in% c("en", "de"))
    .ig_abort("invalid_input", "workbench language must be en or de")
  state$locale <- locale
  state
}
