#' copyToClipboard
#' 
#' Copies then content of a ggraphics widget to the clipboard
#'
#' @param g the ggraphics widget
#'
#' @return
#' 
#' The content of the ggraphics widget on the clipboard, ready to paste
#' 
#' @export
#'
#' @examples
#'  w=gWidgets2::gwindow(visible=FALSE)
#'  gg=gWidgets2::ggraphics(cont=w)
#'  gWidgets2::visible(w) <- TRUE
#'  graphics::hist(rnorm(2000))
#'  copyToClipboard(gg)
#'  # Paste somewhere!
#'  
copyToClipboard = function(g) { 
  da <- getBlock(g)                 # drawing area 
  ww <- da$getAllocation()$allocation$width 
  hh <- da$getAllocation()$allocation$height 
  buf <- gdkPixbufGetFromDrawable(src=da$window, src.x=0, src.y=0,cmap=gdkColormapGetSystem(), 
                                  dest.x=0, dest.y=0, width=ww, height=hh) 
  gtkClipboardGet("CLIPBOARD")$setImage(buf) 
} 
