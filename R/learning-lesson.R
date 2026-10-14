#' Find a bundled worked lesson
#'
#' Returns paths to an original bilingual lesson and its executable synthetic
#' example. Reading this metadata does not execute the example, access the
#' network or write files. The script writes to an explicit output directory
#' when called. External links in the lesson are optional background sources.
#'
#' The Boschloo lesson covers only the two-sided minlike variant for independent
#' binomial groups. It distinguishes numerical optimisation bounds from
#' sampling uncertainty and uses the same variant for prospective planning.
#' A successful teaching example is not full numerical release certification.
#'
#' @param topic Lesson identifier; currently `"boschloo"`.
#' @param locale Language, `"en"` or `"de"`.
#' @return A list with topic, locale, title, exact analysis variant, local
#'   Markdown/HTML/script paths and SHA-256 content hashes. Paths belong to
#'   the installed package and should be treated as read-only.
#' @examples
#' lesson <- learning_lesson("boschloo", "en")
#' basename(lesson$script)
#' @export
learning_lesson <- function(topic = "boschloo", locale = "en") {
  if (!.ig_is_string(topic) || !identical(topic, "boschloo"))
    .ig_abort("invalid_input", "Unknown bundled lesson topic.")
  if (!.ig_is_string(locale) || !locale %in% c("en", "de"))
    .ig_abort("invalid_input", "Lesson locale must be 'en' or 'de'.")
  names <- c(markdown = paste0("boschloo-", locale, ".md"),
             html = paste0("boschloo-", locale, ".html"), script = "boschloo.R")
  paths <- vapply(names, function(name)
    system.file("lessons", "boschloo", name, package = "inferguidR"), "")
  if (any(!nzchar(paths)) || any(!file.exists(paths)))
    .ig_abort("unsupported", "The installed lesson bundle is incomplete.")
  c(list(topic = topic, locale = locale,
    title = if (locale == "en") "Two independent binomial groups: Boschloo's test" else
      "Zwei unabh\u00e4ngige Binomialgruppen: Boschloo-Test",
    analysis_variant_id = "boschloo_two_sided_minlike"), as.list(paths),
    list(sha256 = as.list(vapply(paths, function(path)
      digest::digest(file = path, algo = "sha256"), ""))))
}

.ig_advisor_lesson <- function(recommendation, locale) {
  if (!identical(recommendation$preferred$analysis_variant_id,
                 "boschloo_two_sided_minlike")) return(NULL)
  learning_lesson("boschloo", locale)
}
