#' differentialDeviations
#'
#' Function to see whether deviations differ between groups
#' @param object chromVARDeviations object
#' @param groups either vector of groups or name of column in colData of object 
#' with group information
#' @param alternative only used if there are two groups -- two.sided or one 
#' sided test
#' @param parametric use parametric test. alternatively will use kruskal wallace
#' @return data.frame with p value and adjusted p value
#' @export
#' @author Alicia Schep
#' @importFrom stats aggregate
#' @examples
#' # Load very small example results from computeDeviations
#' data(mini_dev, package = "chromVAR")
#' difdev <- differentialDeviations(mini_dev, "Cell_Type")
differentialDeviations <- function(object, 
                                    groups,
                                    alternative = c("two.sided", "less",
                                                    "greater"), 
                                    parametric = TRUE) {
  stopifnot(is(object,"chromVARDeviations"))
  if (length(groups) == 1 && groups %in% colnames(colData(object))) {
    groups <- colData(object)[[groups]]
  } else if (length(groups) != ncol(object)) {
    stop("invalid groups input, must be vector of lench ncol(object) or column",
         " name from colData(object)")
  }
  
  groups <- as.factor(groups)
  
  alternative <- match.arg(alternative)
  inputs <- deviations(object)
  
  if (parametric) {
    if (nlevels(groups) == 2) {
      # t-test
      p_val <- apply(inputs, 1, t_helper, groups, alternative)
    } else {
      # anova
      p_val <- apply(inputs, 1, anova_helper, groups)
    }
  } else {
    if (nlevels(groups) == 2) {
      # wilcoxon
      p_val <- apply(inputs, 1, wilcoxon_helper, groups, alternative)
    } else {
      # kruskal-wallis
      p_val <- apply(inputs, 1, kw_helper, groups)
    }
  }
  
  p_adj <- p.adjust(p_val, method = "BH")
  return(data.frame(p_value = p_val, p_value_adjusted = p_adj))
}


t_helper <- function(x, groups, alternative) {
  splitx <- split(x, groups)
  return(t.test(splitx[[1]],splitx[[2]],
                alternative = alternative, 
                paired = FALSE,
                var.equal = FALSE)$p.value)
}

anova_helper <- function(x, groups) {
  tmpdf <- data.frame(groups = groups, devs = x)
  res <- oneway.test(devs ~ groups, tmpdf, var.equal = FALSE)
  return(res$p.value)
}

kw_helper <- function(x, groups) {
  tmpdf <- data.frame(groups = groups, devs = x)
  res <- kruskal.test(devs ~ groups, tmpdf)
  return(res$p.value)
}

wilcoxon_helper <- function(x, groups, alternative) {
  splitx <- split(x, groups)
  return(wilcox.test(splitx[[1]], splitx[[2]],
                     alternative = alternative, 
                     paired = FALSE)$p.value)
}

#' differentialVariability
#'
#' Function to determine whether groups differ in variability
#' @param object chromVARDeviations object
#' @param groups either vector of groups or name of column in colData of object 
#' with grouop information
#' @param parametric use parametric test. alternatively will use kruskal wallace
#' @return data.frame with p value and adjusted p value
#' @author Alicia Schep
#' @export
#' @examples
#' # Load very small example results from computeDeviations
#' data(mini_dev, package = "chromVAR")
#' difvar <- differentialVariability(mini_dev, "Cell_Type")
differentialVariability <- function(object, groups, parametric = TRUE) {
  stopifnot(is(object,"chromVARDeviations"))
  if (length(groups) == 1 && groups %in% colnames(colData(object))) {
    groups <- colData(object)[[groups]]
  } else if (length(groups) != ncol(object)) {
    stop("invalid groups input, must be vector of lench ncol(object) or column",
         " name from colData(object)")
  }
  
  groups <- as.factor(groups)
  inputs <- deviationScores(object)
  # Brown-Forsythe test
  if (parametric) {
    p_val <- apply(inputs, 1, bf_var_test, groups)
  } else {
    p_val <- apply(inputs, 1, bf_kw_var_test, groups)
  }
  p_adj <- p.adjust(p_val, method = "BH")
  return(data.frame(p_value = p_val, p_value_adjusted = p_adj))
}

bf_var_test <- function(x, groups) {
  medians <- aggregate(x, list(groups), median)$x
  median_diff <- abs(x - unsplit(medians, groups))
  return(anova(lm(median_diff ~ groups))[1, 5])
}

bf_kw_var_test <- function(x, groups) {
  medians <- aggregate(x, list(groups), median)$x
  median_diff <- abs(x - unsplit(medians, groups))
  return(kruskal.test(median_diff ~ groups)$p.value)
}
