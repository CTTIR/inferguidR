# Strict scalar entry for the workbench; grouping separators are not accepted.
.ig_workbench_parse_number <- function(text, locale = "en") {
  if (!.ig_is_string(locale) || !locale %in% c("en", "de"))
    .ig_abort("invalid_input", "workbench language must be en or de")
  if (!.ig_is_string(text))
    .ig_abort("invalid_input", "enter one number as text")
  text <- trimws(text)
  decimal <- if (locale == "de") "," else "\\."
  pattern <- paste0("^[+-]?(?:[0-9]+(?:", decimal, "[0-9]*)?|", decimal,
                    "[0-9]+)(?:[eE][+-]?[0-9]+)?$")
  if (!grepl(pattern, text, perl = TRUE))
    .ig_abort("invalid_input", if (locale == "de")
      "Eine Zahl mit Dezimalkomma eingeben; keine Tausendertrennzeichen." else
      "Enter one number with a decimal point and no grouping separators.")
  value <- as.double(if (locale == "de") sub(",", ".", text, fixed = TRUE) else text)
  if (!is.finite(value)) .ig_abort("invalid_input", "the entered number must be finite")
  value
}
