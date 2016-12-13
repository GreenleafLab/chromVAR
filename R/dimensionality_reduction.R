

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
  #vars <- vars[ix]

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
      Rtsne::Rtsne(t(mat[which(vars[ix] > input$threshold),]),
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


# deviations_pca <- function(object,
#                             threshold = 1.5,
#                             shiny = interactive(),
#                             ...){
#
#   if (shiny){
#     return(deviations_pca_shiny(object, threshold))
#   } else{
#     vars <- row_sds(assays(object)$z, FALSE)
#     if (threshold > max(vars, na.rm = TRUE)) stop("threshold too high")
#     ix <- which(vars > threshold)
#     mat <- assays(object)[["deviations"]][ix,]
#     ix2 <- remove_correlated_helper(mat, vars[ix])
#     pca_res <- prcomp(t(mat[ix2,]))
#     out <- list(threshold = threshold, pca = pca_res)
#     return(out)
#   }
# }
#
#
# deviations_pca_shiny <- function(object, threshold = 1.5, perplexity = 30){
#
#
#   vars <- row_sds(assays(object)$z, FALSE)
#   mat <- assays(object)$deviations
#
#   ix <- remove_correlated_helper(mat, vars)
#   mat <- mat[ix,]
#   #vars <- vars[ix]
#
#   #if (threshold > max(vars, na.rm = TRUE)) threshold <- quantile(vars[which(vars > 1)], 0.8)
#
#   #starting_value <- median(vars[which(vars > 1)])
#
#   ui <- miniPage(
#     gadgetTitleBar("tsne visualization: adjust parameters on left"),
#     fillRow(miniContentPanel(fillCol(
#       sliderInput("threshold", "Variability threshold", min = 1, max = round(max(vars, na.rm = TRUE),digits = 2),
#                   value = threshold),
#       selectInput("color", "color by", choices = c("none",colnames(colData(object)),
#                                                    rownames(object)[order(vars, decreasing = TRUE)]), selected = "none"))),
#       miniContentPanel(plotlyOutput("plot1", height = "100%")), width = "95%", height ="95%"
#     )
#   )
#
#   server <- function(input, output, session) {
#
#
#     get_pca <- reactive({
#       prcomp(t(mat[which(vars[ix] > input$threshold),]))
#     })
#
#     # Render the plot
#     output$plot1 <- renderPlotly({
#       pca <- get_pca()
#
#       if (input$color == "none"){
#         p1 = ggplot(data.frame(x = tsne$Y[,1], y = tsne$Y[,2], text = colnames(object)),
#                     aes(x = x, y = y,  text = text))  +
#           geom_point(size = 2) + chromVAR_theme(12) +
#           xlab("tSNE dim 1") + ylab("tSNE dim 2")
#       } else if (input$color %in% colnames(colData(object))){
#         p1 = ggplot(data.frame(x = tsne$Y[,1], y = tsne$Y[,2], color = colData(object)[,input$color], text = colnames(object)),
#                     aes(x = x, y = y, col = color, text = text))  +
#           geom_point(size = 2) + chromVAR_theme(12) +
#           labs(col = input$color) +
#           xlab("tSNE dim 1") + ylab("tSNE dim 2") + theme(legend.key.size = grid:::unit(0.5,"lines"))
#       } else{
#         p1 = ggplot(data.frame(x = tsne$Y[,1], y = tsne$Y[,2], color = assays(object[input$color,])[["deviations"]][1,], text = colnames(object)),
#                     aes(x = x, y = y, col = color, text = text)) + geom_point(size= 2) +
#           scale_color_gradient2(name = rowData(object[input$color,])[,"name"],mid = "lightgray", low = "blue", high = "red") + chromVAR_theme(12) +
#           xlab("tSNE dim 1") + ylab("tSNE dim 2") + theme(legend.key.size = grid:::unit(0.5,"lines"))
#       }
#       ggplotly(p1)
#     })
#
#     observeEvent(input$done, {
#       stopApp(list(perplexity = input$perplexity, threshold = input$threshold, tsne = get_tsne()))
#     })
#
#   }
#
#
#   runGadget(ui, server)
#
# }

#' deviations_tsne_motifs
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
deviations_tsne_motifs <- function(motifs = NULL,
                                   kmers = NULL,
                                   threshold = 1.5,
                                   perplexity = 10,
                                   shiny = interactive(),
                                   ...){
  stopifnot(!is.null(motifs) || !is.null(kmers))
  if (shiny){
    return(deviations_tsne_motifs_shiny(motifs = motifs, kmers =kmers, threshold = threshold, perplexity = perplexity))
  } else{
    vars <- row_sds(rbind(assays(motifs)$z,assays(kmers)$z), FALSE)
    if (threshold > max(vars, na.rm = TRUE)) stop("threshold too high")
    ix <- which(vars > threshold)
    mat <- rbind(assays(motifs)$deviations,assays(kmers)$deviations)[ix,]
    tsne_res <- Rtsne::Rtsne(mat, perplexity = perplexity)
    out <- list(threshold = threshold, perplexity = perplexity, tsne = tsne_res, ix = ix)
    return(out)
  }
}


deviations_tsne_motifs_shiny <- function(motifs = NULL, kmers = NULL, threshold = 1.5, perplexity = 30){

  type <- c(rep("motif",nrow(motifs)),rep("kmer",nrow(kmers)))

  vars <- row_sds(rbind(assays(motifs)$z,assays(kmers)$z), FALSE)
  mat <- rbind(assays(motifs)$deviations,assays(kmers)$deviations)

  if (!is.null(motifs)){
    cd <- colData(motifs)
  } else{
    cd <- colData(kmers)
  }

  color_options <- colnames(cd)
  color_options <- color_options[which(sapply(color_options, function(x){
   tmp = cd[,x]
   is.character(tmp) || is.factor(tmp)
  }))]

  ui <- miniPage(
    gadgetTitleBar("tsne visualization: adjust parameters on left"),
    fillRow(miniContentPanel(fillCol(
      sliderInput("perplexity", "Perplexity", min = 3, max = floor(ncol(mat)/2), value = perplexity, round = TRUE),
      sliderInput("threshold", "Variability threshold", min = 1, max = round(max(vars, na.rm = TRUE),digits = 2),
                  value = threshold),
      selectInput("color", "color by", choices = c("none",color_options), selected = "none"))),
      miniContentPanel(plotlyOutput("plot1", height = "100%")), width = "95%", height ="95%"
    )
  )

  server <- function(input, output, session) {


    get_tsne <- reactive({
      Rtsne::Rtsne(mat[which(vars > input$threshold),],
                   perplexity = input$perplexity, check_duplicates = FALSE)
    })

    # Render the plot
    output$plot1 <- renderPlotly({
      tsne <- get_tsne()

      if (input$color == "none"){
        p1 = ggplot(data.frame(x = tsne$Y[,1], y = tsne$Y[,2],
                               text = c(rownames(motifs),rownames(kmers))[which(vars > input$threshold)],
                               Variability = vars[which(vars > input$threshold)],
                               type = type[which(vars > input$threshold)]),
                    aes(x = x, y = y,  text = text, size = Variability, shape = type))  +
          geom_point(size = 2) + chromVAR_theme(12) + scale_size(range=c(0.5,2.5)) + scale_shape_manual(values = c(1,16)) +
          xlab("tSNE dim 1") + ylab("tSNE dim 2")
      } else {
        mean_dev_anno <- sapply(unique(cd[,input$color]),
                                function(x) rowMeans(mat[which(vars > input$threshold),which(cd[,input$color] == x)]))
        motif_anno <- apply(mean_dev_anno,1,function(x) colnames(mean_dev_anno)[which.max(x)])

        p1 = ggplot(data.frame(x = tsne$Y[,1], y = tsne$Y[,2], color = motif_anno,
                               text = c(rownames(motifs),rownames(kmers))[which(vars > input$threshold)],
                               Variability = vars[which(vars > input$threshold)],
                               type = type[which(vars > input$threshold)]),
                    aes(x = x, y = y, col = color, text = text, size = Variability, shape = type))  +
          geom_point(size = 2) + chromVAR_theme(12) + scale_size(range=c(0.5,2.5)) + scale_shape_manual(values = c(1,16)) +
          labs(col = "Max Deviations:\n") +
          xlab("tSNE dim 1") + ylab("tSNE dim 2") + theme(legend.key.size = grid:::unit(0.5,"lines"))
      }
      ggplotly(p1)
    })

    observeEvent(input$done, {
      stopApp(list(perplexity = input$perplexity, threshold = input$threshold, tsne = get_tsne(), ix = which(vars > input$threshold)))
    })

  }


  runGadget(ui, server)

}

