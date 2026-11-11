context("plotting")

# Only plot structure that holds on ggplot2 3.x and 4.x is read: the built
# layer data, the layers, the plot data and the axis labels.
geom_names <- function(p) {
	layers <- p$layers
	skip_if(is.null(layers), "ggplot layers are not accessible")
	vapply(layers, function(l) class(l$geom)[1], "")
}

plot_label <- function(p, aesthetic) {
	if ("get_labs" %in% getNamespaceExports("ggplot2")) {
		ggplot2::get_labs(p)[[aesthetic]]
	} else {
		p$labels[[aesthetic]]
	}
}

y_limits <- function(p) ggplot2::layer_scales(p)$y$get_limits()

# colours drawn for the points, in row order
point_colours <- function(p) ggplot2::ggplot_build(p)$data[[1]]$colour

var_df <- data.frame(name = c("a", "b", "c", "d"),
                     variability = c(1, 4, 2, 3),
                     stringsAsFactors = FALSE)
var_boot <- cbind(var_df,
                  bootstrap_lower_bound = c(0.5, 3, 1, 2.5),
                  bootstrap_upper_bound = c(1.5, 5, 3, 3.5))

test_that("plotVariability without bootstrap bounds plots variability against rank", {
  p <- plotVariability(var_df, n = 0, use_plotly = FALSE)
  expect_is(p, "ggplot")
  expect_equal(geom_names(p), "GeomPoint")
  built <- ggplot2::ggplot_build(p)$data[[1]]
  expect_equal(built$x, c(4, 1, 3, 2))
  expect_equal(built$y, var_df$variability)
  expect_equal(plot_label(p, "x"), "Sorted TFs")
  expect_equal(plot_label(p, "y"), "Variability")
})

test_that("plotVariability y axis starts at zero and extends 5% above the maximum", {
  p <- plotVariability(var_df, n = 0, use_plotly = FALSE)
  expect_equal(y_limits(p), c(0, 4 * 1.05))
  pb <- plotVariability(var_boot, n = 0, use_plotly = FALSE)
  expect_equal(y_limits(pb), c(0, 5 * 1.05))
})

test_that("plotVariability honours the xlab argument", {
  p <- plotVariability(var_df, xlab = "My axis", n = 0, use_plotly = FALSE)
  expect_equal(plot_label(p, "x"), "My axis")
})

test_that("plotVariability with bootstrap bounds adds error bars and maps the bounds", {
  p <- plotVariability(var_boot, n = 0, use_plotly = FALSE)
  expect_equal(geom_names(p), c("GeomPoint", "GeomErrorbar"))
  built <- ggplot2::ggplot_build(p)$data
  expect_equal(built[[2]]$ymin, var_boot$bootstrap_lower_bound)
  expect_equal(built[[2]]$ymax, var_boot$bootstrap_upper_bound)
  expect_equal(built[[1]]$x, c(4, 1, 3, 2))
})

test_that("plotVariability labels the top n ranked annotations", {
  for (n in 1:3) {
    p <- plotVariability(var_boot, n = n, use_plotly = FALSE)
    expect_equal(geom_names(p), c("GeomPoint", "GeomErrorbar", "GeomText"))
    txt <- ggplot2::ggplot_build(p)$data[[3]]
    ranked <- var_boot$name[order(-var_boot$variability)][seq_len(n)]
    expect_equal(sort(as.character(txt$label)), sort(ranked))
    expect_equal(sort(txt$x), seq_len(n))
  }
})

test_that("plotVariability adds no text layer when n is 0", {
  p <- plotVariability(var_boot, n = 0, use_plotly = FALSE)
  expect_false("GeomText" %in% geom_names(p))
})

test_that("plotVariability labels come from the labels argument", {
  p <- plotVariability(var_boot, n = 1, labels = c("w", "x", "y", "z"),
                       use_plotly = FALSE)
  txt <- ggplot2::ggplot_build(p)$data[[3]]
  expect_equal(as.character(txt$label), "x")
})

test_that("plotVariability data keeps the input columns plus rank and annotation", {
  p <- plotVariability(var_boot, n = 0, use_plotly = FALSE)
  expect_equal(p$data$rank, c(4, 1, 3, 2))
  expect_equal(p$data$annotation, var_boot$name)
  expect_equal(p$data$variability, var_boot$variability)
})

test_that("plotVariability ranks tied values with distinct ranks", {
  tied <- data.frame(name = letters[1:3], variability = c(2, 2, 1),
                     stringsAsFactors = FALSE)
  p <- plotVariability(tied, n = 0, use_plotly = FALSE)
  expect_equal(sort(p$data$rank), 1:3)
  expect_equal(p$data$rank[3], 3)
})

# chromVAR_theme -----------------------------------------------------------------

test_that("chromVAR_theme returns a ggplot2 theme with black axis lines and no grid", {
  th <- chromVAR_theme()
  expect_is(th, "theme")
  expect_equal(th$axis.line$colour, "black")
  for (el in c("panel.background", "panel.border", "panel.grid.major",
               "panel.grid.minor", "plot.background", "strip.background",
               "legend.background", "legend.key")) {
    expect_is(th[[el]], "element_blank")
  }
})

test_that("chromVAR_theme passes base_size and base_family on", {
  th <- chromVAR_theme(base_size = 20, base_family = "serif")
  expect_equal(th$text$size, 20)
  expect_equal(th$text$family, "serif")
  th_default <- chromVAR_theme()
  expect_equal(th_default$text$size, 12)
  expect_equal(th_default$text$family, "Helvetica")
})

test_that("chromVAR_theme can be added to a plot", {
  p <- ggplot2::ggplot(data.frame(x = 1:3, y = 1:3), ggplot2::aes(x, y)) +
    ggplot2::geom_point() + chromVAR_theme(8)
  expect_is(p, "ggplot")
  skip_if(is.null(p$theme), "plot theme is not accessible")
  expect_equal(p$theme$text$size, 8)
})

# plotDeviationsTsne -------------------------------------------------------------

make_tsne_obj <- function() {
  dev <- matrix(c(1, 2, 3, 4, 5, 6, 7, 8, 9,
                  9, 8, 7, 6, 5, 4, 3, 2, 1,
                  1, 1, 2, 2, 3, 3, 4, 4, 5,
                  0, 1, 0, 1, 0, 1, 0, 1, 0), nrow = 4, byrow = TRUE,
                dimnames = list(c("m1", "m2", "m3", "m4"), paste0("s", 1:9)))
  z <- dev * 2
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(deviations = dev, z = z),
    colData = S4Vectors::DataFrame(
      type = rep(c("x", "y"), length.out = 9),
      many = paste0("L", 1:9),
      row.names = colnames(dev)),
    rowData = S4Vectors::DataFrame(
      name = c("TF1", "TF2", "TF3", "TF3"), row.names = rownames(dev)))
  new("chromVARDeviations", se)
}

tsne_obj <- make_tsne_obj()
tsne_xy <- cbind(1:9, (9:1) * 10)
rownames(tsne_xy) <- colnames(tsne_obj)

test_that("plotDeviationsTsne returns one plot per sample column, named by column", {
  res <- plotDeviationsTsne(tsne_obj, tsne_xy, sample_column = c("type", "many"),
                            shiny = FALSE)
  expect_equal(names(res), c("type", "many"))
  expect_true(all(vapply(res, ggplot2::is.ggplot, TRUE)))
  d <- res$type$data
  expect_equal(d$x, tsne_xy[, 1], check.attributes = FALSE)
  expect_equal(d$y, tsne_xy[, 2], check.attributes = FALSE)
  expect_equal(as.character(d$color), rep(c("x", "y"), length.out = 9))
  expect_equal(d$text, colnames(tsne_obj))
  expect_equal(geom_names(res$type), "GeomPoint")
  expect_equal(plot_label(res$type, "x"), "tSNE dim 1")
  expect_equal(plot_label(res$type, "y"), "tSNE dim 2")
})

test_that("plotDeviationsTsne uses a named colour scale for up to eight levels only", {
  res <- plotDeviationsTsne(tsne_obj, tsne_xy, sample_column = c("type", "many"),
                            shiny = FALSE)
  dark2 <- RColorBrewer::brewer.pal(8, "Dark2")
  expect_equal(point_colours(res$type), dark2[rep(1:2, length.out = 9)])
  expect_false(any(point_colours(res$many) %in% dark2))
  expect_equal(nlevels(as.factor(res$many$data$color)), 9)

  eight <- tsne_obj
  eight$eight <- rep(letters[1:8], length.out = 9)
  res8 <- plotDeviationsTsne(eight, tsne_xy, sample_column = "eight",
                             shiny = FALSE)
  expect_equal(point_colours(res8$eight), dark2[rep(1:8, length.out = 9)])
})

test_that("plotDeviationsTsne colours annotation plots by deviation z-score", {
  res <- plotDeviationsTsne(tsne_obj, tsne_xy, annotation_name = "m2",
                            shiny = FALSE)
  expect_equal(names(res), "m2")
  expect_equal(res$m2$data$color, unname(SummarizedExperiment::assays(tsne_obj)$z["m2", ]))
  cols <- point_colours(res$m2)
  z <- unname(SummarizedExperiment::assays(tsne_obj)$z["m2", ])
  expect_equal(length(unique(cols)), 9)
  expect_equal(cols[which.max(z)] == cols[which.min(z)], FALSE)
  expect_equal(res$m2$data$x, tsne_xy[, 1], check.attributes = FALSE)
})

test_that("plotDeviationsTsne finds annotations by the name column of rowData", {
  res <- plotDeviationsTsne(tsne_obj, tsne_xy, annotation_name = "TF1",
                            shiny = FALSE)
  expect_equal(names(res), "TF1")
  expect_equal(res$TF1$data$color,
               unname(SummarizedExperiment::assays(tsne_obj)$z["m1", ]))
})

test_that("plotDeviationsTsne picks the highest-variability row among duplicate names", {
  obj <- tsne_obj
  z <- SummarizedExperiment::assays(obj)$z
  z["m3", ] <- c(0, 0, 0, 0, 0, 0, 0, 0, 1)
  z["m4", ] <- rep(c(-5, 5), length.out = 9)
  SummarizedExperiment::assays(obj)$z <- z
  res <- plotDeviationsTsne(obj, tsne_xy, annotation_name = "TF3", shiny = FALSE)
  expect_equal(res$TF3$data$color, unname(z["m4", ]))
})

test_that("plotDeviationsTsne combines sample and annotation plots, sample plots first", {
  res <- plotDeviationsTsne(tsne_obj, tsne_xy, sample_column = "type",
                            annotation_name = c("m1", "TF2"), shiny = FALSE)
  expect_equal(names(res), c("type", "m1", "TF2"))
})

test_that("plotDeviationsTsne accepts the list returned by Rtsne", {
  res <- plotDeviationsTsne(tsne_obj, list(Y = tsne_xy), sample_column = "type",
                            shiny = FALSE)
  expect_equal(res$type$data$x, tsne_xy[, 1], check.attributes = FALSE)
})

test_that("plotDeviationsTsne accepts a data.frame of coordinates", {
  res <- plotDeviationsTsne(tsne_obj, as.data.frame(tsne_xy),
                            sample_column = "type", shiny = FALSE)
  expect_equal(res$type$data$y, tsne_xy[, 2], check.attributes = FALSE)
})

test_that("plotDeviationsTsne returns an empty list when nothing is requested", {
  res <- plotDeviationsTsne(tsne_obj, tsne_xy, shiny = FALSE)
  expect_equal(length(res), 0)
})

test_that("plotDeviationsTsne rejects invalid input", {
  expect_error(plotDeviationsTsne(tsne_obj, tsne_xy[1:4, ],
                                  sample_column = "type", shiny = FALSE),
               "Number of rows of tsne do not match")
  expect_error(plotDeviationsTsne(tsne_obj, tsne_xy,
                                  sample_column = "no_such_column",
                                  shiny = FALSE), "sample_column")
  expect_error(plotDeviationsTsne(tsne_obj, tsne_xy,
                                  annotation_name = "no_such_annotation",
                                  shiny = FALSE), "annotation_name invalid")
})

test_that("plotDeviationsTsne with shiny = TRUE returns an app object without running it", {
  skip_if_not_installed("DT")
  skip_if_not_installed("plotly")
  vdf <- data.frame(variability = c(1, 2, 3, 4), p_value = c(0.1, 0.2, 0.3, 0.4))
  app <- plotDeviationsTsne(tsne_obj, tsne_xy, var_df = vdf,
                            sample_column = "type", shiny = TRUE)
  expect_is(app, "shiny.appobj")
})
