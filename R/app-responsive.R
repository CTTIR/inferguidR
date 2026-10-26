.ig_workbench_plot_labels <- function(plot, width) {
  if (!inherits(plot,"ggplot")) return(plot)
  if (!is.numeric(width) || length(width)!=1L || !is.finite(width) || width<=0) width <- 800
  wrap <- function(text, pixels) {
    if (!is.character(text) || length(text)!=1L) return(text)
    paste(vapply(strsplit(text,"\n",fixed=TRUE)[[1]],function(line)
      paste(strwrap(line,width=max(22L,floor((width-90)/pixels))),collapse="\n"),""),collapse="\n")
  }
  labels <- plot$labels
  plot <- plot + ggplot2::labs(title=wrap(labels$title,12),subtitle=wrap(labels$subtitle,10),
    caption=wrap(labels$caption,8),x=wrap(labels$x,9)) +
    ggplot2::theme(plot.caption=ggplot2::element_text(hjust=0))
  if (width<1000) {
    lattice <- identical(attr(plot,"inferguidR_plot_data")$panel,"lattice")
    plot <- plot + ggplot2::theme(legend.position="bottom",legend.box="vertical") +
      ggplot2::guides(colour=ggplot2::guide_legend(ncol=1),
        linetype=ggplot2::guide_legend(ncol=1),shape=ggplot2::guide_legend(ncol=1),
        size=ggplot2::guide_legend(ncol=1,title.position="top"))
    if(lattice) plot <- plot + ggplot2::guides(fill=ggplot2::guide_colourbar(
      title.position="top",barwidth=grid::unit(4,"cm"))) else
      plot <- plot + ggplot2::guides(fill=ggplot2::guide_legend(ncol=1))
    if(width<500) plot <- plot + ggplot2::theme(plot.margin=ggplot2::margin(8,16,8,8))
  }
  plot
}
