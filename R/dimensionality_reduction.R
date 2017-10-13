

remove_correlated_helper <- function(mat, val, cutoff = 0.9) {
  stopifnot(nrow(mat) == length(val))
  cormat <- cor(t(mat), use = "pairwise.complete.obs")
  diag(cormat) <- NA
  keep <- seq_len(nrow(mat))
  for (i in order(val, decreasing = TRUE)) {
    if (i %in% keep) {
      toremove <- which(cormat[keep, i] >= cutoff)
      if (length(toremove) > 0) 
        keep <- keep[-toremove]
    }
  }
  return(keep)
}


#' deviationsTsne
#'
#' Perform tsne using bias corrected deviations to visualize either cell/sample
#' similarity or motif/kmer/annotation similarity
#' @param object deviations result
#' @param threshold variability threshold -- use only deviations with 
#' variability greater than threshold
#' @param perplexity perplexity parameter for tsne
#' @param theta theta parameter for tsne
#' @param max_iter max iterations parameter for tsne
#' @param what tsne for similarity of samples or annotations?
#' @param shiny load a shiny widget that enable you to explore perplexity and 
#' variability threshold parameter?
#'
#' @return data.frame with two columns for the two dimensions of tSNE output
#' @export
#' @author Alicia Schep
#' @examples 
#' # Load very small example results from computeDeviations
#' data(mini_dev, package = "chromVAR")
#'                          
#' tsne_res <- deviationsTsne(mini_dev, threshold = 0.8, shiny = FALSE)
#' # setting very low variabilitiy threshold because this is mini data set
#' # threshold should generally be above 1
#' # Use plotVariability to get a sense of an appropriate threshold
deviationsTsne <- function(object, 
                            threshold = 1.5, 
                            perplexity = if (what == "samples") 30 else 8, 
                            max_iter = 1000, 
                            theta = 0.5,
                            what = c("samples","annotations"),
                            shiny = FALSE) {
  what <- match.arg(what)
  stopifnot(inherits(object, "chromVARDeviations") || 
              canCoerce(object, "chromVARDeviations"))
  if (what == "samples"){
    if (ncol(object)/3 <= perplexity) {
      message("Perplexity given too high")
      perplexity <- floor((ncol(object) - 1)/3)
      message("Setting perplexity to ", perplexity)
    }
    vars <- row_sds(deviationScores(object), FALSE)
    if (threshold > max(vars, na.rm = TRUE)) 
      stop("Threshold too high")
    if (sum(vars > threshold, na.rm = TRUE) < 2)
      stop("Threshold too high, and/or too few non/NA")
    if (shiny) {
      res <- deviations_tsne_shiny(object, threshold, perplexity, max_iter, 
                                   theta)
      message("Variability threshold used is: ", res$threshold)
      message("Perplexity used is: ", res$perplexity)
      tsne_res <- res$tsne
    } else {
      ix <- which(vars > threshold)
      mat <- deviations(object)[ix, , drop = FALSE]
      ix2 <- remove_correlated_helper(mat, vars[ix])
      tsne_res <- Rtsne::Rtsne(t(mat[ix2, , drop = FALSE]), 
                               perplexity = perplexity, 
                               max_iter = max_iter, 
                               theta = theta)
    }
    out <- tsne_res$Y
    rownames(out) <- colnames(object)
    return(out)
  } else {
    vars <- row_sds(deviationScores(object), FALSE)
    if (threshold > max(vars, na.rm = TRUE)) 
      stop("threshold too high")
    if (sum(vars > threshold, na.rm = TRUE) < 2)
      stop("Threshold too high, and/or too few non/NA")
    if (shiny) {
      res <- deviations_tsne_inv_shiny(object, threshold, perplexity, max_iter, 
                                   theta)
      message("Variability threshold used is: ",res$threshold)
      message("Perplexity used is: ", res$perplexity)
      ix <- which(vars > res$threshold)
      tsne_res <- res$tsne
    } else {
      mat <- deviations(object)
      ix <- which(vars > threshold)
      tsne_res <- Rtsne::Rtsne(mat[ix ,, drop = FALSE], 
                               perplexity = perplexity, 
                               check_duplicates = FALSE,
                               max_iter = max_iter, 
                               theta = theta)
    }
    out <- as.data.frame(tsne_res$Y)
    rownames(out) <- rownames(object)[ix]
    colnames(out) <- c("Dim1","Dim2")
    return(out)  
    }
}

deviations_tsne_shiny <- function(object, threshold = 1.5, perplexity = 30,
                                  max_iter, theta){
  
  
  vars <- row_sds(assays(object)$z, FALSE)
  mat <- assays(object)$deviations
  
  ix <- remove_correlated_helper(mat, vars)
  mat <- mat[ix, , drop = FALSE]
  vars <- vars[ix]

  ui <- miniPage(
    gadgetTitleBar("tsne visualization: adjust parameters on left"),
    fillCol(
      flex = c(1,4),
      miniContentPanel(
        fillRow(
          flex = c(1,1,2),
          numericInput("perplexity", 
                       "Perplexity:", 
                       min = 3, 
                       max = floor(ncol(object)/2), 
                       value = perplexity, 
                       step = 1, width = "90%"),
          numericInput("threshold", 
                       "Variability threshold:", 
                       min = 1, 
                       max = round(max(vars, na.rm = TRUE) - 0.1,
                                   digits = 2),
                       value = threshold, 
                       step = 0.1, width = "90%"),
          selectizeInput("color", 
                         "Color by", 
                         options = list(dropdownParent = 'body'),
                         choices = 
                           c("none",
                             colnames(colData(object)),
                             rownames(object)[order(vars, decreasing = TRUE)]), 
                         selected = "none", width = "90%"))),
      miniContentPanel(plotlyOutput("plot1", height = "100%")), 
      width = "95%", 
      height = "95%")
  )
  
  server <- function(input, output, session) {
    
    
    get_tsne <- reactive({
      Rtsne::Rtsne(t(mat[which(vars > input$threshold), , drop = FALSE]),
                   perplexity = input$perplexity, 
                   max_iter = max_iter, theta = theta)
    })
    
    # Render the plot
    output$plot1 <- renderPlotly({
      tsne <- get_tsne()
      
      if (input$color == "none"){
        p1 <- ggplot(data.frame(x = tsne$Y[,1], y = tsne$Y[,2], 
                               text = colnames(object)),
                    aes_string(x = "x", y = "y",  text = "text"))  +
          geom_point(size = 2) + chromVAR_theme(12) +
          xlab("tSNE dim 1") + ylab("tSNE dim 2")
      } else if (input$color %in% colnames(colData(object))){
        p1 <- ggplot(data.frame(x = tsne$Y[,1], y = tsne$Y[,2], 
                               color = colData(object)[,input$color], 
                               text = colnames(object)),
                    aes_string(x = "x", y = "y", col = "color", 
                               text = "text")) +
          geom_point(size = 2) + chromVAR_theme(12) +
          labs(col = input$color) +
          xlab("tSNE dim 1") + ylab("tSNE dim 2") + 
          theme(legend.key.size = grid::unit(0.5,"lines"))
      } else{
        p1 <- 
          ggplot(data.frame(x = tsne$Y[,1], 
                            y = tsne$Y[,2], 
                            color = deviationScores(object[input$color,])[1,], 
                            text = colnames(object)),
                 aes_string(x = "x", y = "y", col = "color", text = "text")) + 
          geom_point(size= 2) +
          scale_color_gradient2(name = rowData(object[input$color,])[,"name"],
                                mid = "lightgray", low = "blue", high = "red") +
          chromVAR_theme(12) +
          xlab("tSNE dim 1") + ylab("tSNE dim 2") + 
          theme(legend.key.size = grid::unit(0.5,"lines"))
      }
      ggplotly(p1)
    })
    
    observeEvent(input$done, {
      stopApp(list(perplexity = input$perplexity, threshold = input$threshold, 
                   tsne = get_tsne()))
    })
    
  }
  
  
  runGadget(ui, server)
  
}


deviations_tsne_inv_shiny <- function(object, threshold, perplexity,
                                  max_iter, theta){
  
  
  vars <- row_sds(assays(object)$z, FALSE)
  mat <- deviations(object)
  
  ui <- miniPage(
    gadgetTitleBar("tsne visualization: adjust parameters on left"),
    fillCol(flex = c(1,4),miniContentPanel(fillRow(
      numericInput("perplexity", "Perplexity", min = 3, 
                   max = floor(ncol(object)/2), 
                   value = perplexity, step = 1, width = "90%"),
      numericInput("threshold", "Variability threshold:", min = 1, 
                   max = round(max(vars, na.rm = TRUE) - 0.1,
                               digits = 2),
                   value = threshold, step = 0.1, width = "90%"))),
      miniContentPanel(plotlyOutput("plot1", height = "100%")), 
      width = "95%", 
      height ="95%")
  )
  
  server <- function(input, output, session) {
    
    
    get_tsne <- reactive({
      Rtsne::Rtsne(mat[which(vars > input$threshold), , drop = FALSE],
                   perplexity = input$perplexity, 
                   max_iter = max_iter, theta = theta,
                   check_duplicates = FALSE)
    })
    
    # Render the plot
    output$plot1 <- renderPlotly({
      tsne <- get_tsne()
      
      p1 <- 
        ggplot(data.frame(x = tsne$Y[,1], y = tsne$Y[,2], 
                          text = rownames(object)[which(vars > 
                                                          input$threshold)]),
               aes_string(x = "x", y = "y",  text = "text"))  +
        geom_point(size = 2) + chromVAR_theme(12) +
        xlab("tSNE dim 1") + ylab("tSNE dim 2")
      
      ggplotly(p1)
    })
    
    observeEvent(input$done, {
      stopApp(list(perplexity = input$perplexity, threshold = input$threshold, 
                   tsne = get_tsne()))
    })
    
  }
  
  runGadget(ui, server)
  
}

