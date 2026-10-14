# Reduce a finite binary64 argument using 1152 fractional bits of 2/pi.
# Products of 24-bit limbs and the three-limb significand are exact in
# binary64 (each accumulated coefficient is less than 3 * 2^48).
# At the largest finite argument, truncation contributes less than 2^-128
# radians. Small signed remainders avoid subtraction near an integer.
# Native calls are confined to the same small interval used after reduction.
.ig_calc_trig <- function(x, kind) {
  if (!is.finite(x) || abs(x) <= pi / 4) {
    return(switch(kind, sin = sin(x), cos = cos(x), tan = tan(x)))
  }
  base <- 2^24
  inverse <- c(
    10680707, 7228996, 1387004, 2578385, 16069853, 12639074,
    9804092, 4427841, 16666979, 11263675, 12935607, 2387514,
    4345298, 14681673, 3074569, 13734428, 16653803, 1880361,
    10960616, 8533493, 3062596, 8710556, 7349940, 6258241,
    3772886, 3769171, 3798172, 8675211, 12450088, 3874808,
    9961438, 366607, 15675153, 9132554, 7151469, 3571407,
    2607881, 12013382, 4155038, 6285869, 7677882, 13102053,
    15825725, 473591, 9065106, 15363067, 6271263, 9264392
  )
  ax <- abs(x)
  exponent <- min(1023, floor(log2(ax)))
  if (ax / 2^exponent < 1) exponent <- exponent - 1
  significand <- ax * 2^(52 - exponent)
  digits <- numeric(3L)
  for (i in seq_len(3L)) {
    digits[i] <- significand %% base
    significand <- floor(significand / base)
  }
  inverse <- rev(inverse)
  product <- numeric(length(inverse) + 3L)
  for (i in seq_along(digits)) {
    for (j in seq_along(inverse)) {
      product[i + j - 1L] <- product[i + j - 1L] + digits[i] * inverse[j]
    }
  }
  for (i in seq_len(length(product) - 1L)) {
    carry <- floor(product[i] / base)
    product[i] <- product[i] - carry * base
    product[i + 1L] <- product[i + 1L] + carry
  }
  bits <- unlist(lapply(product, function(d) floor(d / 2^(0:23)) %% 2),
                 use.names = FALSE)
  shift <- 1152L + 52L - exponent
  quadrant <- bits[shift + 1L] + 2 * bits[shift + 2L]
  fraction <- bits[seq_len(shift)]
  negative <- fraction[shift] == 1
  if (negative) {
    quadrant <- (quadrant + 1) %% 4
    fraction <- 1 - fraction
    carry <- 1
    for (i in seq_along(fraction)) {
      v <- fraction[i] + carry
      fraction[i] <- v %% 2
      carry <- floor(v / 2)
      if (carry == 0) break
    }
  }
  f <- sum(rev(fraction) * 2^(-seq_len(shift)))
  if (negative) f <- -f
  angle <- f * (pi / 2)
  s <- sin(angle)
  c <- cos(angle)
  sine <- switch(as.integer(quadrant) + 1L, s, c, -s, -c)
  cosine <- switch(as.integer(quadrant) + 1L, c, -s, -c, s)
  if (x < 0) sine <- -sine
  switch(kind, sin = sine, cos = cosine, tan = sine / cosine)
}
