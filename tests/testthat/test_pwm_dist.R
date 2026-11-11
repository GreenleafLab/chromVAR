context("pwmDistance")

library(TFBSTools)

nucs <- c("A", "C", "G", "T")

make_motif <- function(width, a, b) {
  m <- matrix(abs(sin(seq_len(4 * width) * a + b)) + 0.1, nrow = 4,
              dimnames = list(nucs, NULL))
  sweep(m, 2, colSums(m), "/")
}

one_hot <- function(s) {
  letters <- strsplit(s, "")[[1]]
  m <- matrix(0, nrow = 4, ncol = length(letters), dimnames = list(nucs, NULL))
  m[cbind(match(letters, nucs), seq_along(letters))] <- 1
  m
}

rev_comp_matrix <- function(m) {
  out <- m[4:1, rev(seq_len(ncol(m))), drop = FALSE]
  dimnames(out) <- dimnames(m)
  out
}

# Direct implementation of the definition: slide y along x, keep alignments
# overlapping by at least min_overlap columns, average the per-column Euclidean
# distance and scale by sqrt(2) (the largest possible column distance).
ref_dist_one_strand <- function(x, y, min_overlap) {
  n1 <- ncol(x)
  n2 <- ncol(y)
  best <- Inf
  best_offset <- NA
  for (offset in seq(-(n2 - 1), n1 - 1)) {
    lo <- max(0, offset)
    hi <- min(n1, offset + n2) - 1
    if (hi - lo + 1 < min_overlap) next
    cols_x <- lo:hi + 1
    cols_y <- cols_x - offset
    d <- sum(sqrt(colSums((x[, cols_x, drop = FALSE] -
                             y[, cols_y, drop = FALSE])^2))) /
      sqrt(2) / length(cols_x)
    if (d < best) {
      best <- d
      best_offset <- offset
    }
  }
  c(dist = best, offset = best_offset)
}

ref_dist <- function(x, y, min_overlap) {
  fwd <- ref_dist_one_strand(x, y, min_overlap)
  rev <- ref_dist_one_strand(x, rev_comp_matrix(y), min_overlap)
  if (fwd["dist"] <= rev["dist"]) {
    list(dist = fwd[["dist"]], offset = fwd[["offset"]], strand = "+")
  } else {
    list(dist = rev[["dist"]], offset = rev[["offset"]], strand = "-")
  }
}

ref_all <- function(xs, ys, min_overlap) {
  d <- matrix(0, length(xs), length(ys))
  o <- matrix(0, length(xs), length(ys))
  s <- matrix("", length(xs), length(ys))
  for (i in seq_along(xs)) {
    for (j in seq_along(ys)) {
      r <- ref_dist(xs[[i]], ys[[j]], min_overlap)
      d[i, j] <- r$dist
      o[i, j] <- r$offset
      s[i, j] <- r$strand
    }
  }
  list(dist = d, offset = o, strand = s)
}

motifs <- list(make_motif(6, 1.3, 0.2), make_motif(8, 0.9, 1.1),
               make_motif(7, 2.1, 0.5), make_motif(9, 1.7, 2.3))
other_motifs <- list(make_motif(5, 0.6, 0.9), make_motif(10, 2.7, 1.4))

test_that("a single column pair has the hand-computed distance", {
  res <- pwmDistance(list(one_hot("A")), list(one_hot("C")), min_overlap = 1)
  expect_equal(res$dist[1, 1], 1)
  half <- matrix(c(0.5, 0.5, 0, 0), ncol = 1)
  res <- pwmDistance(list(one_hot("A")), list(half), min_overlap = 1)
  expect_equal(res$dist[1, 1], sqrt(0.25 + 0.25) / sqrt(2))
})

test_that("a motif and its reverse complement match with strand '-'", {
  res <- pwmDistance(list(one_hot("A")), list(one_hot("T")), min_overlap = 1)
  expect_equal(res$dist[1, 1], 0)
  expect_equal(res$strand[1, 1], "-")
  res <- pwmDistance(list(one_hot("AC")), list(one_hot("GT")),
                     min_overlap = 2)
  expect_equal(res$dist[1, 1], 0)
  expect_equal(res$strand[1, 1], "-")
  expect_equal(res$offset[1, 1], 0)
})

test_that("a sub-motif is found at the offset where it occurs", {
  res <- pwmDistance(list(one_hot("ACGTTG")), list(one_hot("CGT")),
                     min_overlap = 3)
  expect_equal(res$dist[1, 1], 0)
  expect_equal(res$strand[1, 1], "+")
  expect_equal(res$offset[1, 1], 1)
  res <- pwmDistance(list(one_hot("CGT")), list(one_hot("ACGTTG")),
                     min_overlap = 3)
  expect_equal(res$dist[1, 1], 0)
  expect_equal(res$offset[1, 1], -1)
})

test_that("a sub-motif's reverse complement is found with strand '-'", {
  res <- pwmDistance(list(one_hot("ACGTTG")), list(one_hot("AAC")),
                     min_overlap = 3)
  expect_equal(res$dist[1, 1], 0)
  expect_equal(res$strand[1, 1], "-")
  expect_equal(res$offset[1, 1], 2)
})

test_that("one-list distances match a direct implementation", {
  for (mo in c(3, 5)) {
    res <- pwmDistance(motifs, min_overlap = mo)
    ref <- ref_all(motifs, motifs, mo)
    expect_equal(res$dist, ref$dist)
    expect_equal(res$strand, ref$strand)
    expect_equal(res$offset[lower.tri(ref$offset, diag = TRUE)],
                 ref$offset[lower.tri(ref$offset, diag = TRUE)])
  }
})

test_that("two-list distances match a direct implementation", {
  for (mo in c(2, 5)) {
    res <- pwmDistance(motifs, other_motifs, min_overlap = mo)
    ref <- ref_all(motifs, other_motifs, mo)
    expect_equal(res$dist, ref$dist)
    expect_equal(res$strand, ref$strand)
    expect_equal(res$offset, ref$offset)
  }
})

test_that("distances against reverse-complemented motifs match the reference",
{
  rc_motifs <- lapply(motifs, rev_comp_matrix)
  res <- pwmDistance(motifs, rc_motifs, min_overlap = 4)
  ref <- ref_all(motifs, rc_motifs, 4)
  expect_equal(res$dist, ref$dist)
  expect_equal(res$strand, ref$strand)
  expect_equal(diag(res$dist), rep(0, length(motifs)))
  expect_equal(diag(res$strand), rep("-", length(motifs)))
})

test_that("pwmDistance returns the documented list of matrices", {
  res <- pwmDistance(motifs)
  expect_equal(sort(names(res)), c("dist", "offset", "strand"))
  expect_equal(dim(res$dist), c(4, 4))
  expect_equal(dim(res$offset), c(4, 4))
  expect_equal(dim(res$strand), c(4, 4))
  expect_true(is.numeric(res$dist))
  expect_true(is.integer(res$offset))
  expect_true(is.character(res$strand))
  expect_true(all(res$strand %in% c("+", "-")))
  res2 <- pwmDistance(motifs, other_motifs)
  expect_equal(dim(res2$dist), c(4, 2))
  expect_equal(dim(res2$offset), c(4, 2))
  expect_equal(dim(res2$strand), c(4, 2))
})

test_that("distances lie between 0 and 1", {
  res <- pwmDistance(c(motifs, other_motifs, list(one_hot("AAAAAA"),
                                                  one_hot("CCCCCC"))),
                     min_overlap = 5)
  expect_true(all(res$dist >= 0))
  expect_true(all(res$dist <= 1 + 1e-12))
})

test_that("the distance of a motif to itself is zero at offset 0 on '+'", {
  res <- pwmDistance(motifs, min_overlap = 5)
  expect_equal(diag(res$dist), rep(0, 4))
  expect_equal(diag(res$offset), rep(0, 4))
  expect_equal(diag(res$strand), rep("+", 4))
  res <- pwmDistance(list(one_hot("ACCGTA")), list(one_hot("ACCGTA")),
                     min_overlap = 3)
  expect_equal(res$dist[1, 1], 0)
  expect_equal(res$offset[1, 1], 0)
  expect_equal(res$strand[1, 1], "+")
})

test_that("a motif identical to another has distance 0 to it", {
  res <- pwmDistance(list(motifs[[1]], motifs[[2]], motifs[[1]]))
  expect_equal(res$dist[1, 3], 0)
  expect_equal(res$dist[3, 1], 0)
  expect_equal(res$dist[1, 2], res$dist[3, 2])
})

test_that("the one-list result is symmetric with antisymmetric offsets", {
  res <- pwmDistance(motifs, min_overlap = 4)
  expect_equal(res$dist, t(res$dist))
  expect_equal(res$strand, t(res$strand))
  plus <- res$strand == "+"
  expect_equal(res$offset[plus], -t(res$offset)[plus])
})

test_that("swapping x and y transposes the distances and strands", {
  fwd <- pwmDistance(motifs, other_motifs, min_overlap = 4)
  bwd <- pwmDistance(other_motifs, motifs, min_overlap = 4)
  expect_equal(bwd$dist, t(fwd$dist))
  expect_equal(bwd$strand, t(fwd$strand))
  plus <- fwd$strand == "+"
  expect_equal(bwd$offset[t(plus)], -t(fwd$offset)[t(plus)])
})

test_that("y = x gives the same distances as y = NULL", {
  one <- pwmDistance(motifs, min_overlap = 4)
  two <- pwmDistance(motifs, motifs, min_overlap = 4)
  expect_equal(one$dist, two$dist)
  expect_equal(one$strand, two$strand)
})

test_that("two-list results equal the matching block of the pooled result", {
  pooled <- pwmDistance(c(motifs, other_motifs), min_overlap = 4)
  cross <- pwmDistance(motifs, other_motifs, min_overlap = 4)
  expect_equal(cross$dist, pooled$dist[1:4, 5:6])
  expect_equal(cross$strand, pooled$strand[1:4, 5:6])
})

test_that("distances do not decrease as min_overlap grows", {
  d <- lapply(2:5, function(mo) pwmDistance(motifs, other_motifs,
                                            min_overlap = mo)$dist)
  for (i in 2:length(d)) {
    expect_true(all(d[[i]] >= d[[i - 1]] - 1e-12))
  }
})

test_that("with min_overlap equal to the width, equal-width motifs are aligned", {
  a <- make_motif(6, 1.3, 0.2)
  b <- make_motif(6, 0.7, 1.9)
  res <- pwmDistance(list(a), list(b), min_overlap = 6)
  expect_equal(res$offset[1, 1] , 0)
  direct <- sum(sqrt(colSums((a - b)^2))) / sqrt(2) / 6
  direct_rc <- sum(sqrt(colSums((a - rev_comp_matrix(b))^2))) / sqrt(2) / 6
  expect_equal(res$dist[1, 1], min(direct, direct_rc))
})

test_that("pwmDistance accepts a single matrix", {
  res <- pwmDistance(motifs[[1]], motifs[[2]], min_overlap = 4)
  ref <- ref_dist(motifs[[1]], motifs[[2]], 4)
  expect_equal(res$dist[1, 1], ref$dist)
  expect_equal(res$strand[1, 1], ref$strand)
  expect_equal(res$offset[1, 1], ref$offset)
  expect_equal(dim(pwmDistance(motifs[[1]])$dist), c(1, 1))
})

make_pwm <- function(prob, id, bg = c(A = 0.25, C = 0.25, G = 0.25, T = 0.25),
                     base = 2) {
  mat <- log(prob / bg) / log(base)
  PWMatrix(ID = id, name = id, profileMatrix = mat, bg = bg)
}

test_that("log2-odds PWMatrix input is converted to probabilities", {
  bg <- c(A = 0.3, C = 0.2, G = 0.2, T = 0.3)
  pwms <- PWMatrixList(make_pwm(motifs[[1]], "a", bg),
                       make_pwm(motifs[[2]], "b", bg))
  res <- pwmDistance(pwms, min_overlap = 4)
  ref <- pwmDistance(motifs[1:2], min_overlap = 4)
  expect_equal(res$dist, ref$dist)
  expect_equal(res$strand, ref$strand)
  expect_equal(chromVAR:::pwm_type(pwms[[1]]), "log2")
})

test_that("natural-log-odds PWMatrix input is converted to probabilities", {
  bg <- c(A = 0.3, C = 0.2, G = 0.2, T = 0.3)
  pwm <- make_pwm(motifs[[1]], "a", bg, base = exp(1))
  expect_equal(chromVAR:::pwm_type(pwm), "log")
  res <- pwmDistance(pwm, make_pwm(motifs[[2]], "b", bg, base = exp(1)),
                     min_overlap = 4)
  ref <- pwmDistance(motifs[[1]], motifs[[2]], min_overlap = 4)
  expect_equal(res$dist, ref$dist)
  expect_equal(chromVAR:::pwm_to_prob(pwm)[[1]], motifs[[1]],
               check.attributes = FALSE)
})

test_that("a PWMatrix holding frequencies is used as is", {
  pwm <- PWMatrix(ID = "a", name = "a", profileMatrix = motifs[[1]])
  expect_equal(chromVAR:::pwm_type(pwm), "frequency")
  expect_equal(chromVAR:::pwm_to_prob(pwm)[[1]], motifs[[1]],
               check.attributes = FALSE)
  res <- pwmDistance(PWMatrixList(pwm,
                                  PWMatrix(ID = "b", name = "b",
                                           profileMatrix = motifs[[2]])),
                     min_overlap = 4)
  expect_equal(res$dist, pwmDistance(motifs[1:2], min_overlap = 4)$dist)
})

test_that("PFMatrix counts with unequal column totals are normalised per column", {
	counts <- matrix(c(6, 2, 1, 1, 0, 4, 8, 8), nrow = 4,
	                 dimnames = list(nucs, NULL))
	pf <- PFMatrix(ID = "a", name = "a", profileMatrix = counts)
	expected <- matrix(c(0.6, 0.2, 0.1, 0.1, 0, 0.2, 0.4, 0.4), nrow = 4)
	expect_equal(chromVAR:::pwm_to_prob(pf)[[1]], expected,
	             check.attributes = FALSE)
	expect_equal(chromVAR:::pwm_to_prob(PFMatrixList(pf, pf))[[2]], expected,
	             check.attributes = FALSE)
})

test_that("PFMatrix and PFMatrixList counts are normalised to probabilities", {
  counts1 <- round(motifs[[1]] * 40)
  counts2 <- round(motifs[[2]] * 40)
  counts1[1, ] <- counts1[1, ] + 40 - colSums(counts1)
  counts2[1, ] <- counts2[1, ] + 40 - colSums(counts2)
  pf1 <- PFMatrix(ID = "a", name = "a", profileMatrix = counts1)
  pf2 <- PFMatrix(ID = "b", name = "b", profileMatrix = counts2)
  ref <- pwmDistance(list(counts1 / 40, counts2 / 40), min_overlap = 4)
  expect_equal(pwmDistance(pf1, pf2, min_overlap = 4)$dist[1, 1],
               ref$dist[1, 2])
  res <- pwmDistance(PFMatrixList(pf1, pf2), min_overlap = 4)
  expect_equal(res$dist, ref$dist)
  expect_equal(chromVAR:::pwm_to_prob(pf1)[[1]], counts1 / 40,
               check.attributes = FALSE)
})

test_that("pwm_type recognises frequency, log2 and log matrices", {
  freq <- PWMatrix(ID = "a", name = "a", profileMatrix = motifs[[1]])
  expect_equal(chromVAR:::pwm_type(freq), "frequency")
  expect_equal(chromVAR:::pwm_type(make_pwm(motifs[[1]], "a")), "log2")
  expect_equal(chromVAR:::pwm_type(make_pwm(motifs[[1]], "a", base = exp(1))),
               "log")
})

test_that("a PWMatrix that is not frequencies or log odds is rejected", {
  bad <- PWMatrix(ID = "a", name = "a", profileMatrix = motifs[[1]] * 3)
  expect_error(pwmDistance(bad), "Can't determine format of PWM")
  expect_error(pwmDistance(PWMatrixList(bad)), "Can't determine format")
})

test_that("a plain matrix whose columns do not sum to 1 is rejected", {
  expect_error(pwmDistance(motifs[[1]] * 2), "not equal|TRUE")
  expect_error(pwmDistance(motifs[[1]], motifs[[2]] * 2), "not equal|TRUE")
})

test_that("a list of matrices whose columns do not sum to 1 is rejected", {
  expect_error(pwmDistance(list(motifs[[1]], motifs[[2]] * 2)), "not TRUE")
  expect_error(pwmDistance(motifs[1], list(motifs[[2]] * 2)), "not TRUE")
})

test_that("unsupported input types are rejected", {
  expect_error(pwmDistance(1:3), "incorrect input format")
  expect_error(pwmDistance("ACGT"), "incorrect input format")
  expect_error(pwmDistance(data.frame(a = 1:4)), "incorrect input format")
  expect_error(pwmDistance(motifs, 1:3), "incorrect input format")
})

test_that("distances are the same for permuted input order", {
  perm <- c(3, 1, 4, 2)
  res <- pwmDistance(motifs, min_overlap = 4)
  res_perm <- pwmDistance(motifs[perm], min_overlap = 4)
  expect_equal(res_perm$dist, res$dist[perm, perm])
  expect_equal(res_perm$strand, res$strand[perm, perm])
})

test_that("an empty list yields empty matrices", {
  res <- pwmDistance(list())
  expect_equal(dim(res$dist), c(0, 0))
  res <- pwmDistance(list(), motifs)
  expect_equal(dim(res$dist), c(0, 4))
})

test_that("pwm_euclidean is the normalised sum of column distances", {
  a <- one_hot("AC")
  b <- one_hot("CC")
  expect_equal(chromVAR:::pwm_euclidean(a, b), sqrt(2) / sqrt(2) / 2)
  expect_equal(chromVAR:::pwm_euclidean(a, a), 0)
  expect_equal(chromVAR:::pwm_euclidean(motifs[[1]], motifs[[3]][, 1:6]),
               sum(sqrt(colSums((motifs[[1]] - motifs[[3]][, 1:6])^2))) /
                 sqrt(2) / 6)
})

test_that("pwm_dist_single matches the single-strand reference", {
  res <- chromVAR:::pwm_dist_single(motifs[[2]], motifs[[4]], 4)
  ref <- ref_dist_one_strand(motifs[[2]], motifs[[4]], 4)
  expect_equal(res[1], ref[["dist"]])
  expect_equal(res[2], ref[["offset"]])
})
