
.bl_fmt <- function(x, digits = 6) {
  if (is.null(x) || length(x) == 0L) return("NA")
  ifelse(is.na(x), "NA", formatC(x, digits = digits, format = "g"))
}
