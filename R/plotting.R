# Plotting ---------------------------------------------------------------------

#' plotVariability
#'
#' plot variability of motifs/etc
#' @param variability output from \code{\link{computeVariability}}
#' @param xlab label for x-axis (default is 'Sorted TFs')
#' @param n  number of toppoints to label?
#' @param labels names of sets. if not given, uses rownames of variability
#' @param use_plotly make plot interactive (using plotly)
#' @import ggplot2 
#' @importFrom plotly ggplotly
#' @export
#' @return ggplot or plotly object, depending on whether use_plotly is TRUE
#' @author Alicia Schep
#' @examples
#' # Load very small example results from computeDeviations
#' data(mini_dev, package = "chromVAR")
#' variability <- computeVariability(mini_dev)
#' var_plot <- plotVariability(variability, use_plotly = FALSE) 
plotVariability <- function(variability, xlab = "Sorted TFs", n = 3, 
                             labels = variability$name, 
  use_plotly = interactive()) {

  
  res_df <- cbind(variability, rank = rank(-1 * variability$variability, 
                                           ties.method = "random"), 
    annotation = labels)
  
  ylab <- "Variability"
  
  if ("bootstrap_lower_bound" %ni% colnames(variability)) {
    
    out <- ggplot(res_df, aes_string(x = "rank", y = "variability")) + 
      geom_point() + 
      xlab(xlab) + ylab(ylab) + 
      scale_y_continuous(expand = c(0, 0),
                         limits = 
                           c(0, 
                             max(res_df$variability, na.rm = TRUE) * 1.05)) + 
      chromVAR_theme()
    
  } else {
    
    out <- ggplot(res_df, aes_string(x = "rank", 
                                     y = "variability",
                                     min = "bootstrap_lower_bound", 
                                     max = "bootstrap_upper_bound",
                                     label = "annotation")) + 
      geom_point() + 
      geom_errorbar() + xlab(xlab) + ylab(ylab) + 
      chromVAR_theme() + 
      scale_y_continuous(expand = c(0, 
                                    0), 
                         limits = c(0, 
                                    max(res_df$bootstrap_upper_bound, 
                                        na.rm = TRUE) * 1.05))
  }
  
  if (use_plotly) {
    return(ggplotly(out))
  } else {
    if (n >= 1) {
      top_df <- res_df[res_df$rank <= n, ]
      out <- out + geom_text(data = top_df, size = 3, hjust = -0.45, 
                             col = "Black")
    }
    return(out)
  }
  
}

#' chromVAR_theme
#'
#' theme for use with ggplot2, used by chromVAR plotting functions
#' @param base_size base font size
#' @param base_family base font family
#'
#' @return ggplot2 theme
#' @export 
#' @author Alicia Schep
#' @examples
#' p <- ggplot2::qplot(1:3,1:3) + chromVAR_theme(18)
chromVAR_theme <- function(base_size = 12, base_family="Helvetica"){
  theme_bw(base_size, base_family) %+replace%
    theme(axis.line = element_line(colour = "black"),
          panel.background = element_blank(),
          panel.border = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          plot.background = element_blank(),
          strip.background = element_blank() ,
          legend.background = element_blank(),
          legend.key = element_blank()
    )
}


#' @importFrom DT datatable formatSignif
variability_table <- function(var_df) {
    out <- DT::datatable(var_df, 
                options = list(order = list(2, "desc"), 
                               pageLength = 5), 
    selection = list(mode = "single",
                     selected = which.max(var_df$variability)))  
    formatSignif(out, colnames(var_df), 3)
}


#' plotDeviationsTsne
#'
#' plots sample similarity tsne
#' @param object deviations result object
#' @param tsne result from \code{\link{deviationsTsne}}
#' @param var_df variability result
#' @param sample_column column name for sample data -- colData(object) -- to be
#'  used for coloring points
#' @param annotation_name name of chromVAR annotation for coloring points
#' @param shiny return shiny app?  otherwise return static plots
#' @return shiny app or plots
#' @export
#' @author Alicia Schep
plotDeviationsTsne <- function(object, 
                                 tsne, 
                                 var_df = NULL, 
                                 sample_column = NULL, 
                                 annotation_name = NULL, 
                                 shiny = interactive()) {

  
  if (is.list(tsne) && "Y" %in% names(tsne)) {
    tsne <- tsne$Y
  } 
  
  if (nrow(tsne) != ncol(object)){
    stop("Number of rows of tsne do not match number of columns of object. ", 
         " plotDeviationsTsne takes result of deviationsTsne for samples")
  }
  
  if (shiny) 
    return(plot_deviations_tsne_shiny(object, tsne, var_df, sample_column))
  
  stopifnot(sample_column %in% colnames(colData(object)))
  anno <- colData(object)[, sample_column]
  

    
  
  out <- list()
  if (!is.null(sample_column)) {
    for (i in sample_column) {
      anno <- colData(object)[, i]
      out[[i]] <- ggplot(data.frame(x = tsne[, 1],
                                    y = tsne[, 2], color = anno, 
                                    text = colnames(object)), 
                         aes_string(x = "x", y = "y", 
                                    col = "color", text = "text")) + 
        geom_point(size = 2) + chromVAR_theme() + 
        xlab("tSNE dim 1") + ylab("tSNE dim 2") +
        theme(legend.key.size = grid::unit(0.5, 
                                           "lines"))
      if (nlevels(as.factor(anno)) <= 8){
        out[[i]] <- out[[i]] + scale_color_brewer(palette = "Dark2", name = i)
      } else{
        out[[i]] <- out[[i]] + guides(colour = guide_legend(title = i))
      }
    }
  }
  if (!is.null(annotation_name)) {
    for (i in annotation_name) {
      if (i %in% rownames(object)) {
        ix <- match(c(i), rownames(object))
      } else if (i %in% rowData(object)$name) {
        ix <- which(rowData(object)$name == i)
        if (length(ix) > 1) 
          ix <- ix[which.max(row_sds(assays(object[ix, ])$z))]
      } else if (is.numeric(i)) {
        ix <- i
      } else {
        stop("annotation_name invalid")
      }
      if ("name" %in% colnames(rowData(object))){
        name_val <- rowData(object)[ix, "name"]
      } else{
        name_val <- i
      }
      out[[i]] <- ggplot(data.frame(x = tsne[, 1], y = tsne[, 2], 
                                    color = deviationScores(object)[ix,], 
                                    text = colnames(object)), 
                         aes_string(x = "x", y = "y", col = "color", 
                                    text = "text")) + 
        geom_point(size = 2) + 
        scale_color_gradient2(name = i,
                              mid = "lightgray", low = "blue", high = "red") +
        chromVAR_theme() + 
        xlab("tSNE dim 1") + ylab("tSNE dim 2") + 
        theme(legend.key.size = grid::unit(0.5, 
                                           "lines"))
    }
  }
  return(out)
}

plot_deviations_tsne_shiny <- function(object, tsne, var_df, 
                                       annotation_column) {
  
  if (is.null(var_df)) 
    var_df <- computeVariability(object)
  
  if ("tsne" %in% names(tsne)) 
    tsne <- tsne$tsne
  if (is.null(annotation_column)) {
    annotation_column <- "none"
  } else {
    stopifnot(annotation_column %in% colnames(colData(object)))
  }
  
  
  ui <- fluidPage(
    fluidRow(column(3, h4("Select options for coloring plots:"), 
                    br(), 
                    selectInput("color", "Color first plot by:", 
                                choices = c("none", 
                                            colnames(colData(object))), 
                                selected = annotation_column), 
                    br(), 
                    p(strong("Color second plot by selecting from table:"))), 
             column(9, DT::dataTableOutput("tbl", width = 350))),
    fluidRow(column(6, plotlyOutput("plot1")), 
             column(6, plotlyOutput("plot2"))))
  
  var_tab <- variability_table(var_df)
  
  server <- function(input, output, session) {
    
    
    output$tbl <- DT::renderDataTable(var_tab)
    
    # Render the plot
    output$plot1 <- renderPlotly({
      if (input$color == "none") {
        p1 <- ggplot(data.frame(x = tsne[, 1], 
                                y = tsne[, 2], 
                                text = colnames(object)), 
                     aes_string(x = "x", y = "y", text = "text")) + 
          geom_point(size = 2) + 
          chromVAR_theme(12) + 
          xlab("tSNE dim 1") + ylab("tSNE dim 2") + 
          theme(legend.key.size = grid::unit(0.5, 
                                             "lines"))
      } else {
        p1 <- ggplot(data.frame(x = tsne[, 1], 
                                y = tsne[, 2], 
                                color = colData(object)[,input$color], 
                                text = colnames(object)), 
                     aes_string(x = "x", y = "y", col = "color", 
                                text = "text")) + 
          geom_point(size = 2) + 
          chromVAR_theme(12) + 
          labs(col = input$color) + 
          xlab("tSNE dim 1") + 
          ylab("tSNE dim 2") + 
          theme(legend.key.size = grid::unit(0.5, "lines"))
      }
      ggplotly(p1)
    })
    
    output$plot2 <- renderPlotly({
      s <- input$tbl_rows_selected
      if (length(s) == 0) 
        return(NULL)
      p2 <- 
        ggplot(data.frame(x = tsne[, 1], y = tsne[, 2],
                          color = deviationScores(object)[s,], 
                          text = colnames(object)), 
               aes_string(x = "x", y = "y", col = "color", text = "text")) +
        geom_point(size = 2) + 
        scale_color_gradient2(name = rowData(object)[s,  "name"], 
                              mid = "lightgray", low = "blue", high = "red") +
        chromVAR_theme(12) + 
        xlab("tSNE dim 1") + ylab("tSNE dim 2") + 
        theme(legend.key.size = grid::unit(0.5, "lines"))
      ggplotly(p2)
    })
    
  }
  
  shinyApp(ui, server, options = list(width = 750))
}

