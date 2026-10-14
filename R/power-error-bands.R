# Error regions from retained display rows. Numerical probabilities always
# come from the engine; no area of a sampled curve is used as an estimate.
.ig_power_error_bands <- function(data, critical, ordinate) {
  out <- data[FALSE, , drop = FALSE]
  out$error <- character(); out$group <- integer()
  critical <- unlist(critical, use.names = FALSE)
  critical <- critical[is.finite(critical)]
  group <- 0L
  for (h in c("H0", "H1")) {
    z <- data[data$hypothesis == h, , drop = FALSE]
    z <- z[order(z$x), , drop = FALSE]
    selected <- if (h == "H0") z$reject else !z$reject
    runs <- rle(selected); ends <- cumsum(runs$lengths)
    starts <- ends - runs$lengths + 1L
    for (k in which(runs$values)) {
      idx <- seq.int(starts[k], ends[k])
      if (identical(ordinate, "density")) {
        # A boundary has zero continuous mass. Include its retained density
        # as the polygon edge even when its decision label is rejection.
        neighbours <- c(starts[k] - 1L, ends[k] + 1L)
        neighbours <- neighbours[neighbours >= 1L & neighbours <= nrow(z)]
        idx <- sort(unique(c(idx, neighbours[z$x[neighbours] %in% critical])))
      }
      group <- group + 1L
      rows <- z[idx, , drop = FALSE]
      rows$error <- if (h == "H0") "alpha" else "beta"
      rows$group <- group
      out <- rbind(out, rows)
    }
  }
  rownames(out) <- NULL
  out
}

.ig_power_error_labels <- function(alpha, beta, lang) {
  de <- lang == "de"
  number <- function(x) formatC(x, digits = 5L, format = "g", decimal.mark = if (de) "," else ".")
  c(alpha = paste0(if (de) "Erreichtes alpha = " else "Attained alpha = ", number(alpha)),
    beta = paste0("beta = ", number(beta)))
}

.ig_report_draw_distribution <- function(v, result, locale) {
  d <- do.call(rbind, lapply(v$data, function(row) as.data.frame(row, stringsAsFactors = FALSE)))
  discrete <- identical(v$axes$y, "probability_mass")
  bands <- .ig_power_error_bands(d, v$critical, v$axes$y)
  colours <- c(alpha = "#D98920", beta = "#A9C9C3")
  de <- locale == "de"
  graphics::plot(range(d$x), range(c(0, d$y)), type = "n",
    xlab = if (de) "Teststatistik" else "Test statistic",
    ylab = if (discrete) { if (de) "Wahrscheinlichkeitsmasse" else "Probability mass" } else if (de) "Dichte" else "Density",
    main = if (de) "Verteilungen bei festem Design" else "Fixed-design distributions")
  if (!discrete) for (g in unique(bands$group)) {
    z <- bands[bands$group == g, , drop = FALSE]
    graphics::polygon(c(z$x[1L],z$x,utils::tail(z$x,1L)),c(0,z$y,0),
      col = grDevices::adjustcolor(colours[z$error[1L]], alpha.f = .45), border = NA)
  }
  for (h in c("H0", "H1")) {
    z <- d[d$hypothesis == h, ]; z <- z[order(z$x), ]
    colour <- if (h == "H0") "#304f9c" else "#117777"
    if (discrete) {
      graphics::segments(z$x,0,z$x,z$y,col=colour,lwd=1)
      graphics::points(z$x,z$y,col=colour,pch=if(h=="H0") 1 else 2,cex=.75)
    } else graphics::lines(z$x,z$y,col=colour,lwd=2,lty=if(h=="H0")1 else 2)
  }
  if (discrete && nrow(bands)) graphics::points(bands$x,bands$y,
    pch=ifelse(bands$error=="alpha",22,23),bg=colours[bands$error],col="#444444",cex=1.1)
  critical <- unlist(v$critical,use.names=FALSE)
  critical <- critical[is.finite(critical)]
  if(length(critical)) graphics::abline(v=critical,lty=3,col="#956000")
  graphics::legend("topright",legend=c("H0","H1"),col=c("#304f9c","#117777"),
    lty=if(discrete)c(0,0) else c(1,2),pch=if(discrete)c(1,2) else c(NA,NA),bty="n",cex=.85)
  labels <- .ig_power_error_labels(result$achieved$alpha_attained,result$achieved$beta,locale)
  graphics::legend("topleft",legend=labels,fill=grDevices::adjustcolor(colours,alpha.f=.45),bty="n",cex=.85)
}
