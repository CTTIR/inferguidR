.ig_boschloo_variant_labels <- function(lang) {
  if (!identical(lang, "de")) return(vapply(.bl_variant_info, `[[`, "", "label"))
  c(minlike = "Boschloo: zweiseitig, Minlike (zweiseitiger Fisher-p-Wert)",
    central = "Boschloo: zweiseitig, zentral (doppelter kleinerer einseitiger p-Wert)",
    less = "Boschloo: einseitig, H1: p1 < p2",
    greater = "Boschloo: einseitig, H1: p1 > p2")
}

.ig_power_display_assumptions <- function(result, lang) {
  if (!identical(lang, "de") || !identical(result$method_id, .bl_method_id) ||
      !identical(result$status$code, "ok")) return(result$assumptions)
  variant <- result$spec$options$variant
  ordering <- switch(variant,
    minlike = paste("Ordnung nach dem zweiseitigen Fisher-p-Wert: Summe der bedingten Nullwahrscheinlichkeiten",
      "P(X1 = y | t), die P(X1 = x1 | t) nicht \u00fcbersteigen (relative Bindungstoleranz)."),
    central = paste("Einseitige Fisher-Ordnungen P(X1 <= x1 | t) und P(X1 >= x1 | t);",
      "die beiden unbedingten p-Werte werden getrennt maximiert.",
      "Der zweiseitige p-Wert ist min(1, 2 * min(p_less, p_greater))."),
    less = "Ordnung nach dem einseitigen Fisher-p-Wert P(X1 <= x1 | t); H1: p1 < p2.",
    greater = "Ordnung nach dem einseitigen Fisher-p-Wert P(X1 >= x1 | t); H1: p1 > p2.")
  if (is.null(ordering)) return(result$assumptions)
  c("Zwei unabh\u00e4ngige Binomialstichproben X1 ~ Bin(n1, p1) und X2 ~ Bin(n2, p2) mit festen Gruppengr\u00f6\u00dfen; H0: p1 = p2.",
    paste(.ig_boschloo_variant_labels(lang)[[variant]], ordering),
    paste("Ablehnungsbereich: Tabellen, deren numerische obere Schranke des Boschloo-p-Werts h\u00f6chstens alpha ist.",
      "Testst\u00e4rke: Summe der Binomialwahrscheinlichkeiten in diesem Bereich unter H1.",
      "Erreichtes Niveau: Supremum der Ablehnungswahrscheinlichkeit \u00fcber die gemeinsame Ereigniswahrscheinlichkeit unter H0."),
    "Kein Mid-p, keine Randomisierung, keine Berger-Boos-Einschr\u00e4nkung, keine Nullmarge.")
}
