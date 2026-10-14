# Bracketing, root finding and integer design search ------------------------------
#
# Numerical building blocks of the mode solvers (R/num-modes.R). Every
# routine returns its evidence (brackets, iterations, residuals, the values
# it evaluated) so that results can carry diagnostics, and reports typed
# outcomes instead of plausible numbers when a problem has no answer:
#
# * .ig_scan_roots(): all sign changes of a function on a grid, each refined
#   with uniroot() and classified as a crossing, a steep crossing (continuous,
#   residual limited by the x resolution), an exact grid zero, a jump
#   (discontinuity, e.g. a discrete test) or an unresolved crossing.
# * .ig_bracket_monotone() + .ig_refine_root(): expansion bracketing and
#   refinement for functions that are monotone by construction (power in
#   alpha, the compromise ratio).
# * .ig_integer_search(): smallest admissible integer N whose power reaches
#   a target, with the minimality evidence required by METHOD_CONTRACT
#   (enumeration, or bisection plus predecessor check and a window scan).

# Solver settings (fixed; reported in diagnostics where they matter).
.ig_solver_settings <- list(
  # Largest total sample size searched (integer range of R).
  n_max = .Machine$integer.max - 1,
  # Below this many integers between the smallest admissible design and the
  # solution, minimality is proven by evaluating every admissible design.
  full_scan_span = 4000,
  # Admissible designs below the solution that are re-evaluated after a
  # bisection, to detect non-monotone power.
  window = 50L,
  # Upper limit of a full enumeration (discrete tests, or non-monotone power).
  max_enumeration = 1e6,
  # Largest gap between admissible designs that the stepping accepts.
  max_gap = 1000L,
  # Absolute residual accepted at a continuous root (power units).
  residual_tol = 1e-9,
  # Power differences at or below this are treated as rounding noise of the
  # forward model (the absolute accuracy target of nmath's noncentral
  # routines is 1e-12), for example when deciding whether a sensitivity
  # target equals the minimum of the power curve.
  power_noise = 1e-12,
  # Relative x tolerance of uniroot().
  x_tol = 1e-13
)

# Continuous roots -------------------------------------------------------------------

# Evaluate f safely: errors propagate, non-finite values become NA.
.ig_eval_finite <- function(f, x) {
  v <- f(x)
  if (!is.numeric(v) || length(v) != 1L) stop("objective must return a single number", call. = FALSE)
  if (!is.finite(v)) NA_real_ else as.double(v)
}

# uniroot() on [a, b] (f(a) and f(b) of opposite sign) with a relative x
# tolerance, followed by a residual check. A residual above `residual_tol`
# is classified from the values around the returned root:
#
# * "steep_crossing": f changes sign within +-h of the root and its values
#   grow with the distance from it (|f(x +- 64 h)| at least 8 times
#   |f(x +- h)|). The function is continuous but so steep that the x
#   resolution of double precision cannot bring the residual below the
#   tolerance; the root is accepted and the residual reported.
# * "jump": f changes sign within +-h but keeps the same magnitude further
#   away, i.e. the sign change is a discontinuity (discrete tests).
# * "unresolved": no sign change around the returned value.
.ig_refine_root <- function(f, a, b, fa, fb, residual_tol = .ig_solver_settings$residual_tol,
                            x_tol = .ig_solver_settings$x_tol) {
  tol <- x_tol * max(1, abs(a), abs(b))
  ur <- stats::uniroot(f, c(a, b), f.lower = fa, f.upper = fb, tol = tol, maxiter = 1000L)
  x <- ur$root
  fx <- .ig_eval_finite(f, x)
  kind <- "crossing"
  if (is.na(fx) || abs(fx) > residual_tol) {
    h <- max(10 * tol, 1e-12 * max(1, abs(x)))
    left <- .ig_eval_finite(f, max(a, x - h))
    right <- .ig_eval_finite(f, min(b, x + h))
    sign_change <- !is.na(left) && !is.na(right) && sign(left) != sign(right) &&
      abs(left) > residual_tol && abs(right) > residual_tol
    kind <- "unresolved"
    if (sign_change) {
      far_left <- .ig_eval_finite(f, max(a, x - 64 * h))
      far_right <- .ig_eval_finite(f, min(b, x + 64 * h))
      grows <- !is.na(far_left) && !is.na(far_right) &&
        abs(far_left) >= 8 * abs(left) && abs(far_right) >= 8 * abs(right)
      kind <- if (grows) "steep_crossing" else "jump"
    }
  }
  list(x = x, residual = fx, bracket = c(a, b), f_bracket = c(fa, fb),
       iterations = ur$iter, estim_prec = ur$estim.prec, kind = kind,
       slope = sign(fb - fa))
}

# All roots of f on an increasing grid. Exact zeros at grid points are roots
# of kind "grid_zero"; runs of consecutive zeros are reported as flat regions.
.ig_scan_roots <- function(f, grid, residual_tol = .ig_solver_settings$residual_tol,
                           x_tol = .ig_solver_settings$x_tol) {
  vals <- vapply(grid, function(x) .ig_eval_finite(f, x), numeric(1))
  roots <- list()
  flat <- list()
  n <- length(grid)
  zero <- !is.na(vals) & vals == 0
  i <- 1L
  while (i <= n) {
    if (zero[i]) {
      j <- i
      while (j < n && zero[j + 1L]) j <- j + 1L
      if (j > i) {
        flat[[length(flat) + 1L]] <- c(grid[i], grid[j])
      }
      before <- if (i > 1L) vals[i - 1L] else NA_real_
      after <- if (j < n) vals[j + 1L] else NA_real_
      roots[[length(roots) + 1L]] <- list(
        x = grid[i], residual = 0, bracket = c(grid[i], grid[j]), f_bracket = c(before, after),
        iterations = 0L, estim_prec = 0, kind = "grid_zero",
        slope = if (!is.na(before) && !is.na(after)) sign(after - before) else NA_real_)
      i <- j + 1L
      next
    }
    if (i < n && !is.na(vals[i]) && !is.na(vals[i + 1L]) && !zero[i + 1L] &&
        sign(vals[i]) != sign(vals[i + 1L])) {
      roots[[length(roots) + 1L]] <- .ig_refine_root(f, grid[i], grid[i + 1L], vals[i], vals[i + 1L],
                                                     residual_tol, x_tol)
    }
    i <- i + 1L
  }
  list(roots = roots, flat = flat, grid = grid, values = vals,
       evaluations = length(grid), failed = sum(is.na(vals)))
}

# Bracket the root of a function that is non-increasing (`decreasing = TRUE`)
# or non-decreasing in x, starting at x0 and expanding geometrically within
# (lower, upper). Returns the bracket or a description of why none exists.
.ig_bracket_monotone <- function(f, x0, lower, upper, decreasing, step = 1, max_steps = 400L) {
  sgn <- if (decreasing) -1 else 1 # orient so that g is non-decreasing
  g <- function(x) sgn * .ig_eval_finite(f, x)
  x <- x0
  gx <- g(x)
  trail <- c(x)
  if (is.na(gx)) return(list(found = FALSE, reason = "objective not finite at the start value", trail = trail))
  if (gx == 0) return(list(found = TRUE, exact = TRUE, a = x, b = x, fa = 0, fb = 0, trail = trail))
  dir <- if (gx < 0) 1 else -1
  a <- x
  ga <- gx
  h <- step
  for (k in seq_len(max_steps)) {
    xn <- x + dir * h
    xn <- min(max(xn, lower), upper)
    gn <- g(xn)
    trail <- c(trail, xn)
    if (!is.na(gn) && sign(gn) != sign(ga)) {
      lo <- min(a, xn)
      hi <- max(a, xn)
      flo <- if (lo == a) ga else gn
      fhi <- if (hi == a) ga else gn
      return(list(found = TRUE, exact = FALSE, a = lo, b = hi, fa = sgn * flo, fb = sgn * fhi,
                  trail = trail))
    }
    if (xn == lower || xn == upper) {
      return(list(found = FALSE, reason = "no sign change before the end of the search range",
                  limit = xn, value_at_limit = if (is.na(gn)) NA_real_ else sgn * gn, trail = trail))
    }
    if (!is.na(gn)) {
      a <- xn
      ga <- gn
    }
    x <- xn
    h <- h * 2
  }
  list(found = FALSE, reason = "bracket expansion did not terminate", trail = trail)
}

# Search grids over parameter domains -------------------------------------------------

# A scan grid over a real domain: decades on unbounded sides (0.5-decade
# steps from 1e-6 to 1e4 away from the finite end or from zero), logit
# spacing on bounded intervals. Exclusive bounds are not evaluated.
.ig_domain_grid <- function(domain) {
  lo <- domain$lower
  hi <- domain$upper
  distances <- domain$scan %||% c(1e-6, 1e4)
  dec <- sort(unique(c(distances, 10^seq(log10(distances[1]), log10(distances[2]), by = 0.5))))
  if (is.infinite(lo) && is.infinite(hi)) {
    g <- c(-rev(dec), 0, dec)
  } else if (is.finite(lo) && is.infinite(hi)) {
    g <- c(if (domain$lower_inclusive) lo, lo + dec)
  } else if (is.infinite(lo) && is.finite(hi)) {
    g <- c(hi - rev(dec), if (domain$upper_inclusive) hi)
  } else {
    g <- c(if (domain$lower_inclusive) lo, lo + (hi - lo) * stats::plogis(seq(-20, 20, by = 0.5)),
           if (domain$upper_inclusive) hi)
  }
  g <- sort(unique(g))
  g[g > lo | (domain$lower_inclusive & g == lo)]
}

# Integer design search -----------------------------------------------------------------

# Smallest admissible total sample size N in [n_lo, n_max] whose power is at
# least `target`.
#
# `admissible(N)` says whether the allocation rule yields a design for N;
# `power_at(N)` evaluates the forward model there. Continuous methods
# (`discrete = FALSE`) are searched by doubling and bisection over the
# admissible set; the predecessor design is evaluated and must fail, and
# minimality is proven by enumeration when the solution lies within
# `full_scan_span` integers of the smallest design, otherwise checked on a
# window of designs below it. Any design below the solution that reaches the
# target (non-monotone power) triggers a full enumeration from the smallest
# design. Discrete methods are always enumerated upwards from the smallest
# design (lowest N whose power is not below the target).
.ig_integer_search <- function(power_at, admissible, n_lo, n_max, target, discrete,
                               n_start = NA_real_, settings = .ig_solver_settings) {
  cache <- new.env(parent = emptyenv())
  n_eval <- 0L
  P <- function(N) {
    key <- format(N, scientific = FALSE)
    if (!is.null(cache[[key]])) return(cache[[key]])
    v <- power_at(N)
    if (!is.numeric(v) || length(v) != 1L || is.na(v)) {
      stop(sprintf("the forward model returned no power for N = %s", key), call. = FALSE)
    }
    n_eval <<- n_eval + 1L
    assign(key, v, envir = cache)
    v
  }
  next_up <- function(N) {
    for (k in 0:settings$max_gap) {
      if (N + k > n_max) return(NA_real_)
      if (admissible(N + k)) return(N + k)
    }
    stop(sprintf("no admissible design within %d of N = %s", settings$max_gap,
                 format(N, scientific = FALSE)), call. = FALSE)
  }
  next_down <- function(N, floor_n) {
    for (k in 0:settings$max_gap) {
      if (N - k < floor_n) return(NA_real_)
      if (admissible(N - k)) return(N - k)
    }
    NA_real_
  }
  evaluated <- function(upto = Inf) {
    keys <- ls(cache)
    n <- as.numeric(keys)
    pw <- vapply(keys, function(k) cache[[k]], numeric(1), USE.NAMES = FALSE)
    o <- order(n)
    n <- n[o]
    pw <- pw[o]
    keep <- n <= upto
    list(n = n[keep], power = pw[keep])
  }
  done <- function(N, search, minimality, monotonicity, predecessor = NA_real_, bracket = c(NA, NA),
                   nonmonotone = FALSE, scan = c(NA, NA)) {
    ev <- evaluated(N)
    tail_n <- utils::tail(ev$n, 10L)
    tail_p <- utils::tail(ev$power, 10L)
    diag <- list(
      search = search,
      minimality = minimality,
      monotonicity = monotonicity,
      n_lower = first,
      n_upper_limit = n_max,
      evaluations = n_eval,
      nonmonotone_detected = nonmonotone,
      power_sequence = list(n_total = tail_n, power = tail_p)
    )
    if (!all(is.na(bracket))) diag$bracket <- bracket
    if (!all(is.na(scan))) diag$enumerated_range <- scan
    if (!is.na(predecessor)) {
      diag$predecessor_n_total <- predecessor
      diag$predecessor_power <- P(predecessor)
    }
    list(status = "ok", n = N, power = P(N), diagnostics = diag)
  }
  fail <- function(code, reason, extra = list()) {
    list(status = code, reason = reason,
         diagnostics = c(list(n_lower = first, n_upper_limit = n_max, evaluations = n_eval), extra))
  }
  enumerate <- function(from, to) {
    N <- from
    count <- 0
    while (!is.na(N) && N <= to) {
      if (P(N) >= target) return(N)
      count <- count + 1
      if (count >= settings$max_enumeration) return(NA_real_)
      N <- next_up(N + 1)
    }
    NA_real_
  }

  first <- next_up(n_lo)
  if (is.na(first)) {
    return(fail("infeasible", sprintf("no admissible design with %s <= N <= %s",
                                      format(n_lo, scientific = FALSE), format(n_max, scientific = FALSE))))
  }

  if (discrete) {
    stop_at <- min(n_max, first + settings$max_enumeration)
    N <- enumerate(first, stop_at)
    if (is.na(N)) {
      return(fail("infeasible", sprintf(
        "the target power %s is not reached by any admissible design with N <= %s (enumeration limit)",
        format(target), format(max(evaluated()$n), scientific = FALSE))))
    }
    pred <- next_down(N - 1, first)
    return(done(N, "enumeration", "proven: every admissible design from the smallest one upwards was evaluated",
                "not assumed (discrete test)", predecessor = pred, scan = c(first, N)))
  }

  if (P(first) >= target) {
    return(done(first, "enumeration", "proven: the smallest admissible design reaches the target",
                "not needed"))
  }

  # Doubling: L fails, H reaches the target.
  L <- first
  cand <- if (is.finite(n_start) && n_start > first) ceiling(n_start) else 2 * first
  H <- NA_real_
  repeat {
    C <- next_up(min(cand, n_max))
    if (is.na(C)) C <- next_down(n_max, L + 1)
    if (is.na(C) || C <= L) {
      p_max <- P(L)
      flat <- abs(p_max - P(first)) <= 1e-14
      return(fail("infeasible", sprintf(
        "the target power %s is not reached by any admissible design with N <= %s (power there: %s)%s",
        format(target), format(L, scientific = FALSE), format(signif(p_max, 7)),
        if (flat) "; power does not change with N" else ""),
        list(power_at_limit = p_max, flat = flat)))
    }
    if (P(C) >= target) {
      H <- C
      break
    }
    L <- C
    cand <- max(2 * C, C + 1)
  }

  # Bisection over the admissible set.
  repeat {
    nl <- next_up(L + 1)
    if (is.na(nl) || nl >= H) break
    mid <- floor((L + H) / 2)
    M <- next_up(mid)
    if (is.na(M) || M >= H || M <= L) M <- nl
    if (P(M) >= target) H <- M else L <- M
  }
  bracket <- c(L, H)

  # Minimality evidence.
  if (H - first <= settings$full_scan_span) {
    N <- enumerate(first, H)
    nonmono <- !is.na(N) && N < H
    return(done(N, "bisection", "proven: every admissible design below the solution was evaluated",
                if (nonmono) "violated: a smaller design reaches the target (non-monotone power)" else
                  "consistent with monotone power on the enumerated range",
                predecessor = next_down(N - 1, first), bracket = bracket, nonmonotone = nonmono,
                scan = c(first, N)))
  }
  window_ok <- TRUE
  W <- L
  for (k in seq_len(settings$window)) {
    W <- next_down(W - 1, first)
    if (is.na(W)) break
    if (P(W) >= target) {
      window_ok <- FALSE
      break
    }
  }
  if (window_ok) {
    return(done(H, "bisection",
                sprintf("predecessor fails; the %d admissible designs below it also fail", settings$window),
                "assumed (continuous method) and consistent with the window scan",
                predecessor = L, bracket = bracket))
  }
  if (H - first > settings$max_enumeration) {
    return(fail("numerical_failure", sprintf(
      "power is not monotone in N (N = %s reaches the target below the bisection result) and the range %s..%s is too large to enumerate",
      format(W, scientific = FALSE), format(first, scientific = FALSE), format(H, scientific = FALSE)),
      list(bracket = bracket)))
  }
  N <- enumerate(first, H)
  done(N, "enumeration", "proven: every admissible design below the solution was evaluated",
       "violated: power is not monotone in N", predecessor = next_down(N - 1, first),
       bracket = bracket, nonmonotone = TRUE, scan = c(first, N))
}
