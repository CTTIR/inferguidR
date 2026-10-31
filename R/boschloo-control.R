
# Shared input checks ---------------------------------------------------------------

.bl_is_count <- function(x, lower = 0) {
  is.numeric(x) && length(x) == 1L && is.finite(x) && x == round(x) && x >= lower
}

.bl_check_alpha <- function(alpha, name = "alpha") {
  if (!is.numeric(alpha) || length(alpha) != 1L || !is.finite(alpha) || alpha <= 0 || alpha >= 1) {
    .ig_abort_invalid(sprintf("`%s` must be a number in (0, 1)", name))
  }
  as.double(alpha)
}

.bl_check_prob <- function(p, name) {
  if (!is.numeric(p) || length(p) != 1L || !is.finite(p) || p < 0 || p > 1) {
    .ig_abort_invalid(sprintf("`%s` must be a probability in [0, 1]", name))
  }
  as.double(p)
}
