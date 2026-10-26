# Copyright (c) 2026 Raban Heller. MIT licence; see ../../../LICENSE.md.
# Independent restricted-likelihood score bisection and full-table enumeration.
ref_pmf <- function(n, p) exp(lchoose(n, 0:n) + (0:n) * log(p) + (n - 0:n) * log1p(-p))

ref_Phibar <- function(z) {
    half <- 0.5 * stats::pgamma(z^2/2, 0.5, lower.tail = FALSE)
    ifelse(z >= 0, half, 1 - half)
}

ref_Phi <- function(z) ref_Phibar(-z)

ref_bisect <- function(f, lo, hi, n = 200) {
    flo <- f(lo)
    for (i in seq_len(n)) {
        m <- (lo + hi)/2
        if (m == lo || m == hi)
            break
        fm <- f(m)
        if (sign(fm) == sign(flo)) {
            lo <- m
            flo <- fm
        }
        else hi <- m
    }
    (lo + hi)/2
}

ref_qnorm_upper <- function(a) ref_bisect(function(z) ref_Phibar(z) - a, -40, 40)

ref_t_density <- function(x, df) {
    exp(lgamma((df + 1)/2) - lgamma(df/2) - 0.5 * log(df * pi) - (df + 1)/2 * log1p(x^2/df))
}

ref_t_upper <- function(t, df) {
    if (t < 0)
        return(1 - ref_t_upper(-t, df))
    stats::integrate(ref_t_density, t, Inf, df = df, rel.tol = 1e-12, subdivisions = 1000L)$value
}

ref_qt_upper <- function(a, df) ref_bisect(function(t) ref_t_upper(t, df) - a, -60, 60, n = 120)

ref_chisq_density <- function(v, df) exp((df/2 - 1) * log(v) - v/2 - (df/2) * log(2) - lgamma(df/2))

ref_v_upper <- function(df) df + 60 * sqrt(2 * df) + 200

ref_nct_upper <- function(c, df, ncp) {
    stats::integrate(function(v) ref_Phibar(c * sqrt(v/df) - ncp) * ref_chisq_density(v, df), 0, ref_v_upper(df),
        rel.tol = 1e-11, subdivisions = 2000L)$value
}

ref_nct_lower <- function(c, df, ncp) {
    stats::integrate(function(v) ref_Phi(c * sqrt(v/df) - ncp) * ref_chisq_density(v, df), 0, ref_v_upper(df),
        rel.tol = 1e-11, subdivisions = 2000L)$value
}

ref_term <- function(cnt, prob) ifelse(cnt == 0, 0, cnt/prob)

ref_rmle <- function(x1, x2, n1, n2, null) {
    sc <- null$scale
    h <- null$value
    p1_of <- switch(sc, difference = function(p) p + h, ratio = function(p) h * p, odds_ratio = function(p) h *
        p/(1 - p + h * p))
    dp1 <- switch(sc, difference = function(p) rep(1, length(p)), ratio = function(p) rep(h, length(p)),
        odds_ratio = function(p) h/(1 - p + h * p)^2)
    lo0 <- switch(sc, difference = max(0, -h), ratio = 0, odds_ratio = 0)
    hi0 <- switch(sc, difference = min(1, 1 - h), ratio = min(1, 1/h), odds_ratio = 1)
    score <- function(p) {
        a <- p1_of(p)
        (ref_term(x1, a) - ref_term(n1 - x1, 1 - a)) * dp1(p) + ref_term(x2, p) - ref_term(n2 - x2, 1 -
            p)
    }
    eps <- 1e-15
    lo <- rep(lo0, length(x1))
    hi <- rep(hi0, length(x1))
    s_lo <- score(pmin(lo + eps, hi0))
    s_hi <- score(pmax(hi - eps, lo0))
    at_lo <- s_lo <= 0
    at_hi <- s_hi >= 0
    for (i in 1:200) {
        m <- (lo + hi)/2
        s <- score(m)
        up <- s > 0
        lo[up] <- m[up]
        hi[!up] <- m[!up]
    }
    p2 <- (lo + hi)/2
    p2[at_lo] <- lo0
    p2[at_hi & !at_lo] <- hi0
    cbind(p1_of(p2), p2)
}

ref_T <- function(stat, x1, x2, n1, n2, null) {
    N <- n1 + n2
    m <- x1 + x2
    cross <- x1 * n2 - x2 * n1
    fin <- function(num, den) {
        out <- num/den
        out[num == 0] <- 0
        i <- num != 0 & den == 0
        out[i] <- sign(num[i]) * Inf
        out
    }
    ss <- x1 * (n1 - x1)/n1 + x2 * (n2 - x2)/n2
    if (null$scale == "difference" && null$value == 0 && !startsWith(stat, "score_")) {
        pooled_den <- sqrt(n1 * n2 * m * (N - m)/N)
        unp_den <- sqrt(n2^2 * x1 * (n1 - x1)/n1 + n1^2 * x2 * (n2 - x2)/n2)
        cc <- sign(cross) * pmax(abs(cross) - N/2, 0)
        return(switch(stat, z_pooled = fin(cross, pooled_den), z_pooled_cc = fin(cc, pooled_den), z_unpooled = fin(cross,
            unp_den), z_unpooled_cc = fin(cc, unp_den), mantel_haenszel = sign(cross) * sqrt(ifelse(m ==
            0 | m == N, 0, (N - 1) * cross^2/(n1 * n2 * m * (N - m)))), likelihood_ratio = {
            e <- function(o, ex) ifelse(o == 0, 0, o * log(o/ex))
            g2 <- 2 * (e(x1, n1 * m/N) + e(n1 - x1, n1 * (N - m)/N) + e(x2, n2 * m/N) + e(n2 - x2, n2 *
                (N - m)/N))
            sign(cross) * sqrt(pmax(g2, 0))
        }, t_test = fin(cross/(n1 * n2), sqrt(ss/(N - 2) * (1/n1 + 1/n2))), stop("no inequality statistic ",
            stat)))
    }
    if (null$scale == "difference") {
        d0 <- null$value
        num <- x1/n1 - x2/n2 - d0
        if (stat %in% c("z_pooled", "z_pooled_cc", "z_unpooled", "z_unpooled_cc", "t_test")) {
            pb <- m/N
            ccn <- sign(num) * pmax(abs(num) - (1/n1 + 1/n2)/2, 0)
            unp <- sqrt(x1 * (n1 - x1)/n1^3 + x2 * (n2 - x2)/n2^3)
            return(switch(stat, z_pooled = fin(num, sqrt(pb * (1 - pb) * (1/n1 + 1/n2))), z_pooled_cc = fin(ccn,
                sqrt(pb * (1 - pb) * (1/n1 + 1/n2))), z_unpooled = fin(num, unp), z_unpooled_cc = fin(ccn,
                unp), t_test = fin(num, sqrt(ss/(N - 2) * (1/n1 + 1/n2)))))
        }
        phi <- 1
    }
    else if (null$scale == "ratio") {
        num <- x1/n1 - null$value * x2/n2
        phi <- null$value
    }
    else {
        r <- ref_rmle(x1, x2, n1, n2, null)
        num <- x1 - n1 * r[, 1]
        num[abs(num) < 1e-09] <- 0
        v <- 1/(1/(n1 * r[, 1] * (1 - r[, 1])) + 1/(n2 * r[, 2] * (1 - r[, 2])))
        v[!is.finite(v)] <- 0
        if (stat == "score_miettinen_nurminen")
            v <- v * N/(N - 1)
        return(fin(num, sqrt(v)))
    }
    r <- ref_rmle(x1, x2, n1, n2, null)
    a <- r[, 1]
    b <- r[, 2]
    v <- a * (1 - a)/n1 + phi^2 * b * (1 - b)/n2
    z <- fin(num, sqrt(if (stat == "score_miettinen_nurminen") v * N/(N - 1) else v))
    if (stat != "score_gart_nam")
        return(z)
    g <- ifelse(v > 0, (a * (1 - a) * (1 - 2 * a)/n1^2 - phi^3 * b * (1 - b) * (1 - 2 * b)/n2^2)/v^1.5,
        0)
    S <- vapply(seq_along(z), function(i) {
        zi <- z[i]
        gi <- g[i]
        if (!is.finite(zi))
            return(zi)
        if (v[i] == 0)
            return(0)
        if (abs(gi) < 1e-300)
            return(zi)
        f <- function(s) s + gi * (s^2 - 1)/6 - zi
        turn <- -3/gi
        if (gi > 0) {
            if (f(turn) > 0)
                return(NaN)
            hi <- max(turn + 1, zi + 1)
            while (f(hi) < 0) hi <- hi * 2 + 1
            ref_bisect(f, turn, hi)
        }
        else {
            if (f(turn) < 0)
                return(NaN)
            lo <- min(turn - 1, zi - 1)
            while (f(lo) > 0) lo <- lo * 2 - 1
            ref_bisect(f, lo, turn)
        }
    }, numeric(1))
    S
}

ref_crit <- function(stat, alpha, tails, N) {
    a <- if (tails == "two")
        alpha/2
    else alpha
    if (stat == "t_test")
        ref_qt_upper(a, N - 2)
    else ref_qnorm_upper(a)
}

ref_exact <- function(stat, n1, n2, p1, p10, p2, null, alpha, tails, dir) {
    g <- expand.grid(x1 = 0:n1, x2 = 0:n2)
    Tv <- ref_T(stat, g$x1, g$x2, n1, n2, null)
    cv <- ref_crit(stat, alpha, tails, n1 + n2)
    rej <- if (tails == "two")
        (Tv >= cv | Tv <= -cv)
    else if (dir > 0)
        Tv >= cv
    else Tv <= -cv
    rej[is.na(rej)] <- FALSE
    P1 <- as.vector(outer(ref_pmf(n1, p1), ref_pmf(n2, p2)))
    P0 <- as.vector(outer(ref_pmf(n1, p10), ref_pmf(n2, p2)))
    c(power = sum(P1[rej]), beta = sum(P1[!rej]), alpha = sum(P0[rej]))
}

ref_large <- function(stat, n1, n2, p1, p2, null, alpha, tails, dir) {
    N <- n1 + n2
    k <- 1/n1 + 1/n2
    two <- tails == "two"
    c <- ref_crit(stat, alpha, tails, N)
    sideA <- function(d, s0, s1, cc = 0) {
        pw <- ref_Phi((dir * d - cc - c * s0)/s1)
        if (two)
            pw <- pw + ref_Phi((-dir * d - cc - c * s0)/s1)
        pw
    }
    sideB <- function(d, s0) {
        pw <- ref_Phi(dir * d/s0 - c)
        if (two)
            pw <- pw + ref_Phi(-dir * d/s0 - c)
        pw
    }
    sc <- null$scale
    h <- null$value
    if (sc == "difference") {
        d <- p1 - p2 - h
        s1 <- sqrt(p1 * (1 - p1)/n1 + p2 * (1 - p2)/n2)
        pb <- (n1 * p1 + n2 * p2)/N
        sp <- sqrt(pb * (1 - pb) * k)
        if (stat %in% c("z_pooled", "likelihood_ratio"))
            return(sideA(d, sp, s1))
        if (stat == "z_pooled_cc")
            return(sideA(d, sp, s1, k/2))
        if (stat == "z_unpooled")
            return(sideA(d, s1, s1))
        if (stat == "z_unpooled_cc")
            return(sideA(d, s1, s1, k/2))
        if (stat == "mantel_haenszel")
            return(sideB(d, sp * sqrt(N/(N - 1))))
        if (stat == "t_test") {
            st <- sqrt((n1 * p1 * (1 - p1) + n2 * p2 * (1 - p2))/(N - 2) * k)
            ncp <- dir * d/st
            pw <- ref_nct_upper(c, N - 2, ncp)
            if (two)
                pw <- pw + ref_nct_lower(-c, N - 2, ncp)
            return(pw)
        }
    }
    r <- ref_rmle(n1 * p1, n2 * p2, n1, n2, null)
    mn <- if (stat == "score_miettinen_nurminen")
        sqrt(N/(N - 1))
    else 1
    if (sc == "odds_ratio") {
        v <- 1/(1/(n1 * r[1, 1] * (1 - r[1, 1])) + 1/(n2 * r[1, 2] * (1 - r[1, 2])))
        return(sideB(n1 * (p1 - r[1, 1]), sqrt(v) * mn))
    }
    phi <- if (sc == "ratio")
        h
    else 1
    d <- if (sc == "ratio")
        p1 - h * p2
    else p1 - p2 - h
    s1 <- sqrt(p1 * (1 - p1)/n1 + phi^2 * p2 * (1 - p2)/n2)
    if (stat == "score_gart_nam") {
        S <- ref_T("score_gart_nam", n1 * p1, n2 * p2, n1, n2, null)
        return(sideA(d, d/S, s1))
    }
    s0 <- sqrt(r[1, 1] * (1 - r[1, 1])/n1 + phi^2 * r[1, 2] * (1 - r[1, 2])/n2) * mn
    sideA(d, s0, s1)
}
