context("differential tests")

make_dev <- function(dev, z = dev, group = NULL) {
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(deviations = dev, z = z),
    colData = S4Vectors::DataFrame(
      grp = if (is.null(group)) rep("a", ncol(dev)) else group,
      row.names = colnames(dev)))
  new("chromVARDeviations", se)
}

local_seed(42)
ref_dev <- matrix(rnorm(3 * 12), nrow = 3,
                  dimnames = list(c("r1", "r2", "r3"), paste0("s", 1:12)))
ref_z <- matrix(rnorm(3 * 12), nrow = 3, dimnames = dimnames(ref_dev))
grp2 <- rep(c("A", "B"), each = 6)
grp3 <- rep(c("A", "B", "C"), each = 4)
obj2 <- make_dev(ref_dev, ref_z, grp2)
obj3 <- make_dev(ref_dev, ref_z, grp3)

row_apply <- function(mat, f) vapply(seq_len(nrow(mat)), function(i) f(mat[i, ]), 1)

test_that("differentialDeviations with two groups matches Welch t-tests on deviations", {
  for (alt in c("two.sided", "less", "greater")) {
    res <- differentialDeviations(obj2, "grp", alternative = alt)
    expected <- row_apply(ref_dev, function(x)
      t.test(x[grp2 == "A"], x[grp2 == "B"], alternative = alt)$p.value)
    expect_equal(res$p_value, expected)
    expect_equal(res$p_value_adjusted, p.adjust(expected, method = "BH"))
  }
})

test_that("differentialDeviations output is a data.frame keyed by annotation", {
  res <- differentialDeviations(obj2, "grp")
  expect_is(res, "data.frame")
  expect_equal(colnames(res), c("p_value", "p_value_adjusted"))
  expect_equal(rownames(res), rownames(ref_dev))
  expect_true(all(res$p_value >= 0 & res$p_value <= 1))
  expect_true(all(res$p_value_adjusted >= res$p_value - 1e-12))
})

test_that("differentialDeviations alternatives are consistent with each other", {
  two <- differentialDeviations(obj2, "grp", "two.sided")$p_value
  less <- differentialDeviations(obj2, "grp", "less")$p_value
  greater <- differentialDeviations(obj2, "grp", "greater")$p_value
  expect_equal(less + greater, rep(1, 3))
  expect_equal(two, 2 * pmin(less, greater))
})

test_that("differentialDeviations uses the first factor level as the reference group", {
  swapped <- ifelse(grp2 == "A", "B", "A")
  less_swapped <- differentialDeviations(obj2, swapped, "less")$p_value
  greater <- differentialDeviations(obj2, "grp", "greater")$p_value
  expect_equal(less_swapped, greater)
})

test_that("differentialDeviations accepts a vector, a column name, or a factor", {
  by_col <- differentialDeviations(obj2, "grp")
  expect_equal(differentialDeviations(obj2, grp2), by_col)
  expect_equal(differentialDeviations(obj2, factor(grp2)), by_col)
  expect_equal(differentialDeviations(obj2, rep(c(1, 2), each = 6)), by_col)
})

test_that("differentialDeviations gives hand-computed Welch p-value for one separated case", {
  dev <- matrix(c(1, 2, 3, 4, 5, 7), nrow = 1,
                dimnames = list("x", paste0("s", 1:6)))
  obj <- make_dev(dev, group = rep(c("g1", "g2"), each = 3))
  # group means 2 and 16/3, sample variances 1 and 7/3
  se2 <- 1 / 3 + (7 / 3) / 3
  tstat <- (2 - 16 / 3) / sqrt(se2)
  df <- se2^2 / ((1 / 3)^2 / 2 + ((7 / 3) / 3)^2 / 2)
  expect_equal(differentialDeviations(obj, "grp")$p_value,
               2 * pt(-abs(tstat), df))
})

test_that("differentialDeviations ignores z and uses the deviations assay", {
  obj_a <- make_dev(ref_dev, ref_z, grp2)
  obj_b <- make_dev(ref_dev, matrix(0, 3, 12, dimnames = dimnames(ref_dev)),
                    grp2)
  expect_equal(differentialDeviations(obj_a, "grp"),
               differentialDeviations(obj_b, "grp"))
})

test_that("differentialDeviations with three groups matches Welch one-way ANOVA", {
  res <- differentialDeviations(obj3, "grp")
  expected <- row_apply(ref_dev, function(x)
    oneway.test(x ~ factor(grp3), var.equal = FALSE)$p.value)
  expect_equal(res$p_value, expected)
  expect_equal(res$p_value_adjusted, p.adjust(expected, method = "BH"))
})

test_that("differentialDeviations alternative is ignored with more than two groups", {
  expect_equal(differentialDeviations(obj3, "grp", alternative = "less"),
               differentialDeviations(obj3, "grp", alternative = "greater"))
})

test_that("differentialDeviations non-parametric two-group matches wilcox.test", {
  for (alt in c("two.sided", "less", "greater")) {
    res <- differentialDeviations(obj2, "grp", alternative = alt,
                                  parametric = FALSE)
    expected <- row_apply(ref_dev, function(x)
      wilcox.test(x[grp2 == "A"], x[grp2 == "B"], alternative = alt)$p.value)
    expect_equal(res$p_value, expected)
    expect_equal(res$p_value_adjusted, p.adjust(expected, method = "BH"))
  }
})

test_that("differentialDeviations Wilcoxon p-values for fully separated groups are exact", {
  dev <- matrix(c(1, 2, 3, 4, 5, 6), nrow = 1,
                dimnames = list("x", paste0("s", 1:6)))
  obj <- make_dev(dev, group = rep(c("g1", "g2"), each = 3))
  # W = 0 of 9 pairs; 20 equally likely rank assignments
  expect_equal(differentialDeviations(obj, "grp", "two.sided",
                                      parametric = FALSE)$p_value, 2 / 20)
  expect_equal(differentialDeviations(obj, "grp", "less",
                                      parametric = FALSE)$p_value, 1 / 20)
  expect_equal(differentialDeviations(obj, "grp", "greater",
                                      parametric = FALSE)$p_value, 1)
})

test_that("differentialDeviations non-parametric three-group matches kruskal.test", {
  res <- differentialDeviations(obj3, "grp", parametric = FALSE)
  expected <- row_apply(ref_dev, function(x)
    kruskal.test(x ~ factor(grp3))$p.value)
  expect_equal(res$p_value, expected)
})

test_that("differentialDeviations Kruskal-Wallis p-value for separated groups is hand-computed", {
  dev <- matrix(1:9, nrow = 1, dimnames = list("x", paste0("s", 1:9)))
  obj <- make_dev(dev, group = rep(c("g1", "g2", "g3"), each = 3))
  # rank sums 6, 15, 24: H = 12/(9*10) * (36 + 225 + 576) / 3 - 30 = 7.2,
  # chi-square with 2 df has survival function exp(-H / 2)
  expect_equal(differentialDeviations(obj, "grp", parametric = FALSE)$p_value,
               exp(-3.6))
})

test_that("differentialDeviations works on the mini_dev example", {
  data(mini_dev, package = "chromVAR")
  res <- differentialDeviations(mini_dev, "Cell_Type")
  devs <- SummarizedExperiment::assays(mini_dev)$deviations
  g <- mini_dev$Cell_Type
  expected <- row_apply(devs, function(x)
    t.test(x[g == "GM"], x[g == "H1"])$p.value)
  expect_equal(res$p_value, expected)
  expect_equal(rownames(res), rownames(mini_dev))
})

test_that("differentialDeviations rejects invalid input", {
  expect_error(differentialDeviations(ref_dev, grp2), "is\\(object")
  expect_error(differentialDeviations(obj2, grp2[1:5]),
               "invalid groups input")
  expect_error(differentialDeviations(obj2, "not_a_column"),
               "invalid groups input")
  expect_error(differentialDeviations(obj2, "grp", alternative = "bogus"),
               "should be one of")
})

# differentialVariability -------------------------------------------------------

bf_reference <- function(x, g, test = c("anova", "kruskal")) {
  g <- factor(g)
  med <- tapply(x, g, median)
  absdev <- abs(x - med[as.integer(g)])
  if (match.arg(test) == "anova") {
    oneway.test(absdev ~ g, var.equal = TRUE)$p.value
  } else {
    kruskal.test(absdev ~ g)$p.value
  }
}

test_that("differentialVariability matches Brown-Forsythe computed from definition", {
  for (obj_g in list(list(obj2, grp2), list(obj3, grp3))) {
    res <- differentialVariability(obj_g[[1]], "grp")
    expected <- row_apply(ref_z, function(x) bf_reference(x, obj_g[[2]]))
    expect_equal(res$p_value, expected)
    expect_equal(res$p_value_adjusted, p.adjust(expected, method = "BH"))
  }
})

test_that("differentialVariability hand-computed two-group case", {
  z <- matrix(c(1, 2, 3, 0, 5, 10), nrow = 1,
              dimnames = list("x", paste0("s", 1:6)))
  obj <- make_dev(dev = matrix(0, 1, 6, dimnames = dimnames(z)), z = z,
                  group = rep(c("g1", "g2"), each = 3))
  # absolute deviations from medians: (1, 0, 1) and (5, 0, 5)
  # SSB = 32 / 3, SSW = 52 / 3, F = (32 / 3) / ((52 / 3) / 4)
  fstat <- (32 / 3) / ((52 / 3) / 4)
  expect_equal(differentialVariability(obj, "grp")$p_value,
               pf(fstat, 1, 4, lower.tail = FALSE))
})

test_that("differentialVariability non-parametric matches Kruskal-Wallis on absolute deviations", {
  res <- differentialVariability(obj3, "grp", parametric = FALSE)
  expected <- row_apply(ref_z, function(x) bf_reference(x, grp3, "kruskal"))
  expect_equal(res$p_value, expected)
  res2 <- differentialVariability(obj2, "grp", parametric = FALSE)
  expect_equal(res2$p_value,
               row_apply(ref_z, function(x) bf_reference(x, grp2, "kruskal")))
})

test_that("differentialVariability uses z scores, not deviations, and returns a keyed data.frame", {
  obj_b <- make_dev(matrix(0, 3, 12, dimnames = dimnames(ref_dev)), ref_z, grp2)
  res <- differentialVariability(obj_b, grp2)
  expect_equal(res, differentialVariability(obj2, "grp"))
  expect_is(res, "data.frame")
  expect_equal(colnames(res), c("p_value", "p_value_adjusted"))
  expect_equal(rownames(res), rownames(ref_z))
})

test_that("differentialVariability is invariant to shifting z within the whole sample", {
  shifted <- make_dev(ref_dev, ref_z + 100, grp2)
  expect_equal(differentialVariability(shifted, "grp")$p_value,
               differentialVariability(obj2, "grp")$p_value)
})

test_that("differentialVariability detects a variance difference", {
  local_seed(1)
  z <- rbind(c(rnorm(10, sd = 0.1), rnorm(10, sd = 10)),
             rnorm(20))
  dimnames(z) <- list(c("wide", "same"), paste0("s", 1:20))
  obj <- make_dev(z, z, rep(c("lo", "hi"), each = 10))
  res <- differentialVariability(obj, "grp")
  expect_lt(res["wide", "p_value"], 0.01)
  expect_equal(res["wide", "p_value"],
               bf_reference(z["wide", ], rep(c("lo", "hi"), each = 10)))
})

test_that("differentialVariability rejects invalid input", {
  expect_error(differentialVariability(ref_z, grp2), "is\\(object")
  expect_error(differentialVariability(obj2, grp2[1:5]),
               "invalid groups input")
  expect_error(differentialVariability(obj2, "missing_column"),
               "invalid groups input")
})
