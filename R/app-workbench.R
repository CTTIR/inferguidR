#' Create the interactive study workbench
#'
#' Creates the bundled Shiny application without starting a server. The
#' current research preview links a registry-driven planning workspace, synthetic
#' classifier examples and the guided method advisor. Input controls follow the
#' method definitions, including option applicability and solve targets.
#' Pending numerical domains are not presented as
#' released or externally validated calculations.
#'
#' @return A Shiny application object. Requires 'shiny' and 'ggplot2'.
#' @seealso [launch_app()]
#' @export
workbench_app <- function() {
  for (package in c("shiny", "ggplot2")) if (!requireNamespace(package, quietly = TRUE))
    .ig_abort("unsupported", sprintf("The workbench requires the suggested package '%s'.", package))
  shiny::addResourcePath("inferguidr-assets", system.file("www", package = "inferguidR"))
  ui <- shiny::fluidPage(
    shiny::tags$head(shiny::tags$style(shiny::HTML(paste(
      "body{font-family:'DejaVu Sans',Arial,sans-serif;color:#152330;background:#fff;font-size:16px}",
      ".ig-header{padding:22px 12px 16px;border-bottom:2px solid #334E9A;display:flex;justify-content:space-between;gap:24px}",
      ".ig-header h1{font-size:28px;margin:0 0 8px}.ig-header p{max-width:75ch;margin:0}",
      ".ig-study{display:grid;grid-template-columns:minmax(210px,1fr) minmax(360px,3fr) minmax(210px,1fr);gap:26px;padding:22px 10px}",
      ".ig-controls{border-right:1px solid #DCE4EA;padding-right:20px}.ig-inspector{border-left:1px solid #DCE4EA;padding-left:20px;overflow-wrap:anywhere}",
      "input,select,button{font-size:16px!important}.btn{border-radius:3px}.btn:focus,a:focus,input:focus,select:focus{outline:3px solid #956000!important;outline-offset:2px}",
      ".ig-recommendation{max-width:85ch;margin:24px 12px}.ig-recommendation h4{margin-top:24px;font-size:19px}",
      "fieldset{max-width:78ch;margin:22px 0}legend{font-size:22px}details{margin:16px 0}table{font-variant-numeric:tabular-nums}",
      "@media(max-width:1100px){.ig-study{grid-template-columns:240px 1fr}.ig-inspector{grid-column:1/-1;border-left:0;border-top:1px solid #DCE4EA;padding:16px 0}}",
      "@media(max-width:720px){.ig-study{display:block}.ig-controls{border-right:0;border-bottom:1px solid #DCE4EA;padding:0 0 20px}.ig-header{display:block}.ig-figure{padding-top:16px}}"
    )))),
    shiny::tags$header(class = "ig-header", shiny::tags$div(shiny::tags$h1(
      shiny::tags$img(src = "inferguidr-assets/logo.svg", alt = "",
                      height = "58", style = "margin-right:12px;vertical-align:middle"),
      "inferguidR"),
      shiny::tags$p("Study design, evidence and decisions / Studiendesign, Evidenz und Entscheidungen")),
      shiny::selectInput("language", "Language / Sprache", c("English" = "en", "Deutsch" = "de"), width = "180px")),
    shiny::tabsetPanel(id = "workspace",
      shiny::tabPanel("Plan a study / Studie planen", value = "power", .ig_power_workspace_ui("power")),
      shiny::tabPanel("Evaluate a classifier / Klassifikator", value = "classifier", .ig_classifier_workspace_ui("classifier")),
      shiny::tabPanel("Choose a method / Methode finden", value = "advisor", .ig_advisor_ui("advisor"))))
  server <- function(input, output, session) {
    language <- shiny::reactive(input$language %||% "en")
    power <- .ig_power_workspace_server("power", language)
    classifier <- .ig_classifier_workspace_server("classifier", language)
    advisor <- .ig_advisor_server("advisor", language)
  }
  shiny::shinyApp(ui, server)
}
