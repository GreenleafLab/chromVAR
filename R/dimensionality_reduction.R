

remove_correlated_helper <- function(mat, val, cutoff = 0.9){
  stopifnot(nrow(mat) == length(val))
  cormat <- cor(t(mat), use = "pairwise.complete.obs")
  diag(cormat) <- NA
  keep <- 1:nrow(mat)
  for (i in order(val, decreasing = TRUE)){
    if (i %in% keep){
      toremove <- which(cormat[keep,i] >= cutoff)
      if (length(toremove) > 0 ) keep = keep[-toremove]
    }
  }
  return(keep)
}


#' deviations_tsne
#'
#' @param object deviations result
#' @param threshold variability threshold -- use only deviations with variability greater than threshold
#' @param perplexity perplexity parameter for tsne
#' @param shiny output shiny gadget for choosing parameters & visualizing results
#' @param ... 
#'
#' @return list with three elements: threshold used, perplexity used, tsne results 
#' from \code{\link[Rtsne]{Rtsne}}
#' @export
#'
#' @examples
deviations_tsne <- function(object,
                            threshold = 1.5, 
                            perplexity = 30,
                            shiny = interactive(),
                            ...){
  
  if (shiny){
    return(deviations_tsne_shiny(object, threshold, perplexity))
  } else{
    vars <- row_sds(assays(object)$z, FALSE)
    if (threshold > max(vars, na.rm = TRUE)) stop("threshold too high")
    ix <- which(vars > threshold)
    mat <- assays(object)[["deviations"]][ix,]
    ix2 <- remove_correlated_helper(mat, vars[ix])
    tsne_res <- Rtsne::Rtsne(t(mat[ix2,]), perplexity = perplexity)
    out <- list(threshold = threshold, perplexity = perplexity, tsne = tsne_res)
    return(out)
  }
}


deviations_tsne_shiny <- function(object, threshold = 1.5, perplexity = 30){
  
  
  vars <- row_sds(assays(object)$z, FALSE)
  mat <- assays(object)$deviations
  
  ix <- remove_correlated_helper(mat, vars)
  mat <- mat[ix,]
  vars <- vars[ix]
  
  #if (threshold > max(vars, na.rm = TRUE)) threshold <- quantile(vars[which(vars > 1)], 0.8)
  
  #starting_value <- median(vars[which(vars > 1)])
  
  ui <- miniPage(
    gadgetTitleBar("tsne visualization: adjust parameters on left"),
    fillRow(miniContentPanel(fillCol(
      sliderInput("perplexity", "Perplexity", min = 3, max = floor(ncol(object)/2), value = perplexity, round = TRUE),
      sliderInput("threshold", "Variability threshold", min = 1, max = round(max(vars, na.rm = TRUE),digits = 2), 
                  value = threshold),
      selectInput("color", "color by", choices = c("none",colnames(colData(object)),
                                                   rownames(object)[order(vars, decreasing = TRUE)]), selected = "none"))),
      miniContentPanel(plotlyOutput("plot1", height = "100%")), width = "95%", height ="95%"
    )
  )
  
  server <- function(input, output, session) {
    
    
    get_tsne <- reactive({
      Rtsne::Rtsne(t(mat[which(vars > input$threshold),]),
                   perplexity = input$perplexity)
    })
    
    # Render the plot
    output$plot1 <- renderPlotly({
      tsne <- get_tsne()
      
      if (input$color == "none"){
        p1 = ggplot(data.frame(x = tsne$Y[,1], y = tsne$Y[,2], text = colnames(object)),
                    aes(x = x, y = y,  text = text))  + 
          geom_point(size = 2) + chromVAR_theme(12) + 
          xlab("tSNE dim 1") + ylab("tSNE dim 2")
      } else if (input$color %in% colnames(colData(object))){
        p1 = ggplot(data.frame(x = tsne$Y[,1], y = tsne$Y[,2], color = colData(object)[,input$color], text = colnames(object)),
                    aes(x = x, y = y, col = color, text = text))  + 
          geom_point(size = 2) + chromVAR_theme(12) + 
          labs(col = input$color) + 
          xlab("tSNE dim 1") + ylab("tSNE dim 2") + theme(legend.key.size = grid:::unit(0.5,"lines"))
      } else{
        p1 = ggplot(data.frame(x = tsne$Y[,1], y = tsne$Y[,2], color = assays(object[input$color,])[["deviations"]][1,], text = colnames(object)),
                    aes(x = x, y = y, col = color, text = text)) + geom_point(size= 2) +
          scale_color_gradient2(name = rowData(object[input$color,])[,"name"],mid = "lightgray", low = "blue", high = "red") + chromVAR_theme(12) +
          xlab("tSNE dim 1") + ylab("tSNE dim 2") + theme(legend.key.size = grid:::unit(0.5,"lines"))
      }
      ggplotly(p1)
    })
    
    observeEvent(input$done, {
      stopApp(list(perplexity = input$perplexity, threshold = input$threshold, tsne = get_tsne()))
    })
    
  }
  
  
  runGadget(ui, server)
  
}
