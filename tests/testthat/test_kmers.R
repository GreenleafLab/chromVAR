context("kmers")

library(Biostrings)
library(SummarizedExperiment)
library(TFBSTools)

nucs <- c("A", "C", "G", "T")

rc_string <- function(x) {
  vapply(strsplit(chartr("ACGT", "TGCA", x), ""),
         function(y) paste(rev(y), collapse = ""), "")
}

all_kmers <- function(k) {
  apply(expand.grid(rep(list(nucs), k),
                    stringsAsFactors = FALSE)[, rev(seq_len(k)), drop = FALSE],
        1, paste, collapse = "")
}

canon <- function(x) pmin(x, rc_string(x))

canonical_kmers <- function(k) sort(unique(canon(all_kmers(k))))

hamming <- function(a, b) {
  sum(strsplit(a, "")[[1]] != strsplit(b, "")[[1]])
}

# All k-mers sharing a k - s long suffix/prefix with the seed, for a shift of
# s positions to the left (-s) or right (+s).
brute_shifted <- function(kmer, s, side) {
  k <- nchar(kmer)
  km <- all_kmers(k)
  if (side == "left") {
    km[substr(km, s + 1, k) == substr(kmer, 1, k - s)]
  } else {
    km[substr(km, 1, k - s) == substr(kmer, s + 1, k)]
  }
}

brute_overlaps <- function(kmer, max_extend, dir = "both") {
  out <- character(0)
  for (s in seq_len(max_extend)) {
    if (dir %in% c("both", "left")) out <- c(out, brute_shifted(kmer, s, "left"))
    if (dir %in% c("both", "right")) out <- c(out, brute_shifted(kmer, s, "right"))
  }
  out
}

make_cov_mat <- function(k) {
  kmers <- canonical_kmers(k)
  n <- length(kmers)
  m <- matrix(sin(seq_len(n * n) * 0.37) * 0.1, n, n,
              dimnames = list(kmers, kmers))
  diag(m) <- 1
  m
}

cov6 <- make_cov_mat(6)

test_that("get_reverse_complement reverse complements strings", {
  expect_equal(chromVAR:::get_reverse_complement(c("AAC", "GGT", "ACGT")),
               c("GTT", "ACC", "ACGT"))
  expect_equal(chromVAR:::get_reverse_complement(rc_string(canonical_kmers(4))),
               canonical_kmers(4)[order(canonical_kmers(4))])
})

test_that("seq_to_pwm encodes nucleotides as indicator columns", {
  m <- chromVAR:::seq_to_pwm("ACGTA")
  expect_equal(dim(m), c(4, 5))
  expect_equal(rownames(m), nucs)
  expected <- matrix(0, 4, 5, dimnames = list(nucs, NULL))
  expected[cbind(c(1, 2, 3, 4, 1), 1:5)] <- 1
  expect_equal(m, expected)
})

test_that("seq_to_pwm expands IUPAC codes and applies the mismatch value", {
  m <- chromVAR:::seq_to_pwm("RN")
  expect_equal(unname(m[, 1]), c(1, 0, 1, 0))
  expect_equal(unname(m[, 2]), c(1, 1, 1, 1))
  m <- chromVAR:::seq_to_pwm("CY", mismatch = 0.1)
  expect_equal(unname(m[, 1]), c(0.1, 1, 0.1, 0.1))
  expect_equal(unname(m[, 2]), c(0.1, 1, 0.1, 1))
  expect_equal(unname(chromVAR:::seq_to_pwm("W")[, 1]), c(1, 0, 0, 1))
})

test_that("get_mm_kmers lists all single-nucleotide mismatches", {
  res <- chromVAR:::get_mm_kmers("ACG")
  expect_equal(length(res), 9)
  expect_equal(length(unique(res)), 9)
  expect_true(all(vapply(res, hamming, 0, b = "ACG") == 1))
  expect_equal(res, c("CCG", "GCG", "TCG", "AAG", "AGG", "ATG", "ACA", "ACC",
                      "ACT"))
  res <- chromVAR:::get_mm_kmers("TTTTT")
  expect_equal(length(res), 15)
})

test_that("expected_n_overlap_kmers sums 4^i over the extensions", {
  expect_equal(vapply(0:4, chromVAR:::expected_n_overlap_kmers, 0),
               c(0, 4, 20, 84, 340))
})

test_that("get_overlap_kmers matches a brute-force definition", {
  for (kmer in c("ACG", "ACGT", "TTAGCA")) {
    for (ext in 1:2) {
      res <- chromVAR:::get_overlap_kmers(kmer, ext)
      expect_equal(sort(res), sort(brute_overlaps(kmer, ext)))
      expect_equal(length(res), 2 * chromVAR:::expected_n_overlap_kmers(ext))
    }
  }
  expect_equal(sort(chromVAR:::get_overlap_kmers("ACGTA", 3)),
               sort(brute_overlaps("ACGTA", 3)))
})

test_that("get_overlap_kmers respects the direction", {
  for (ext in 1:3) {
    left <- chromVAR:::get_overlap_kmers("ACGTA", ext, "left")
    right <- chromVAR:::get_overlap_kmers("ACGTA", ext, "right")
    expect_equal(sort(left), sort(brute_overlaps("ACGTA", ext, "left")))
    expect_equal(sort(right), sort(brute_overlaps("ACGTA", ext, "right")))
    expect_equal(length(left), chromVAR:::expected_n_overlap_kmers(ext))
  }
  expect_equal(chromVAR:::get_overlap_kmers("ACG", 1, "left"),
               c("AAC", "CAC", "GAC", "TAC"))
  expect_equal(chromVAR:::get_overlap_kmers("ACG", 1, "right"),
               c("CGA", "CGC", "CGG", "CGT"))
})

test_that("get_overlap_kmers returns nothing for no extension", {
  expect_equal(length(chromVAR:::get_overlap_kmers("ACG", 0)), 0)
  expect_equal(length(chromVAR:::get_overlap_kmers("ACG", -1)), 0)
})

test_that("kmers_to_names maps kmers to the orientation present in names", {
  expect_equal(chromVAR:::kmers_to_names(c("AAA", "TTT", "ACG", "CGT"),
                                         c("AAA", "ACG")),
               c("AAA", "AAA", "ACG", "ACG"))
})

test_that("deviationsCovariability divides covariances by the row variance", {
  z <- matrix(sin(seq_len(6 * 12) * 0.9) + cos(seq_len(6 * 12) * 0.31), 6, 12,
              dimnames = list(paste0("k", 1:6), NULL))
  res <- deviationsCovariability(SummarizedExperiment(assays = list(z = z)))
  expect_equal(dim(res), c(6, 6))
  expect_equal(dimnames(res), list(paste0("k", 1:6), paste0("k", 1:6)))
  ref <- matrix(0, 6, 6)
  for (i in 1:6) {
    for (j in 1:6) {
      zi <- z[i, ]
      zj <- z[j, ]
      ref[i, j] <- sum((zi - mean(zi)) * (zj - mean(zj))) /
        sum((zi - mean(zi))^2)
    }
  }
  expect_equal(unname(res), ref)
  expect_equal(unname(diag(res)), rep(1, 6))
})

test_that("deviationsCovariability of scaled rows is asymmetric in the scale", {
  base <- sin(seq_len(10) * 1.7)
  other <- cos(seq_len(10) * 0.6)
  z <- rbind(base, 2 * base + other)
  res <- deviationsCovariability(SummarizedExperiment(assays = list(z = z)))
  expect_equal(res[1, 2], 2 + sum((base - mean(base)) * (other - mean(other))) /
                 sum((base - mean(base))^2))
  expect_equal(res[1, 1], 1)
  expect_equal(res[2, 2], 1)
})

neighborhood <- function(seed, ext) {
  mm <- unlist(lapply(seq_len(nchar(seed)), function(i) {
    vapply(setdiff(nucs, substr(seed, i, i)), function(n) {
      paste0(substr(seed, 1, i - 1), n, substr(seed, i + 1, nchar(seed)))
    }, "")
  }), use.names = FALSE)
  list(mm = mm, o = brute_overlaps(seed, ext))
}

test_that("get_covariable_kmers lists the seed, mismatches and overlaps", {
  kg <- chromVAR:::get_covariable_kmers("ACGTAC", cov6)
  nb <- neighborhood("ACGTAC", 2)
  expect_equal(colnames(kg), c("kmer", "mismatch", "shift", "covariability",
                               "pval", "pval.adj"))
  expect_equal(nrow(kg), 1 + 18 + 40)
  expect_equal(kg$kmer[1], "ACGTAC")
  expect_equal(sort(kg$kmer[2:19]), sort(nb$mm))
  expect_equal(sort(kg$kmer[20:59]), sort(nb$o))
  expect_true(is.na(kg$mismatch[1]))
  expect_equal(kg$mismatch[2:19], rep(1:6, each = 3))
  expect_equal(kg$shift[1:19], rep(0, 19))
  expect_equal(kg$shift[20:59], c(rep(-1, 4), rep(1, 4), rep(-2, 16),
                                  rep(2, 16)))
})

test_that("each row of get_covariable_kmers has the stated mismatch or shift", {
  kg <- chromVAR:::get_covariable_kmers("ACGTAC", cov6)
  for (i in which(!is.na(kg$mismatch))) {
    pos <- which(strsplit(kg$kmer[i], "")[[1]] !=
                   strsplit("ACGTAC", "")[[1]])
    expect_equal(pos, kg$mismatch[i])
  }
  for (i in which(kg$shift < 0)) {
    s <- -kg$shift[i]
    expect_equal(substr(kg$kmer[i], s + 1, 6), substr("ACGTAC", 1, 6 - s))
  }
  for (i in which(kg$shift > 0)) {
    s <- kg$shift[i]
    expect_equal(substr(kg$kmer[i], 1, 6 - s), substr("ACGTAC", s + 1, 6))
  }
})

test_that("get_covariable_kmers covariability is read from the matrix", {
  kg <- chromVAR:::get_covariable_kmers("ACGTAC", cov6)
  expect_equal(kg$covariability, unname(cov6["ACGTAC", canon(kg$kmer)]))
  expect_equal(kg$covariability[1], 1)
})

test_that("get_covariable_kmers p values follow the null distribution", {
  kg <- chromVAR:::get_covariable_kmers("ACGTAC", cov6)
  nb <- neighborhood("ACGTAC", 2)
  omm <- unique(unlist(lapply(nb$mm, brute_overlaps, max_extend = 3)))
  excluded <- unique(canon(c("ACGTAC", nb$mm, nb$o, omm)))
  null_vals <- cov6["ACGTAC", setdiff(colnames(cov6), excluded)]
  cands <- unique(canon(c(nb$mm, nb$o)))
  scores <- (cov6["ACGTAC", cands] - mean(null_vals)) / sd(null_vals)
  pvals <- pnorm(scores, lower.tail = FALSE)
  padj <- p.adjust(pvals)
  expect_equal(kg$pval, unname(pvals[canon(kg$kmer)]))
  expect_equal(kg$pval.adj, unname(padj[canon(kg$kmer)]))
  expect_true(all(kg$pval.adj >= kg$pval - 1e-12))
  expect_true(all(kg$pval.adj <= 1))
})

test_that("get_covariable_kmers canonicalises the seed", {
  kg <- chromVAR:::get_covariable_kmers("ACGTAC", cov6)
  rc_seed <- rc_string("ACGTAC")
  expect_false(rc_seed %in% colnames(cov6))
  kg_rc <- chromVAR:::get_covariable_kmers(rc_seed, cov6)
  expect_equal(kg_rc, kg)
})

test_that("get_covariable_kmers max_extend sets the overlap rows", {
  kg0 <- chromVAR:::get_covariable_kmers("ACGTAC", cov6, max_extend = 0)
  expect_equal(nrow(kg0), 1 + 18)
  expect_equal(kg0$shift, rep(0, 19))
  kg1 <- chromVAR:::get_covariable_kmers("ACGTAC", cov6, max_extend = 1)
  expect_equal(nrow(kg1), 1 + 18 + 8)
  expect_equal(kg1$shift[20:27], c(rep(-1, 4), rep(1, 4)))
  expect_equal(sort(kg1$kmer[20:27]),
               sort(brute_overlaps("ACGTAC", 1)))
})

test_that("get_null_kmer_dist excludes the seed's neighbourhood", {
  nd <- chromVAR:::get_null_kmer_dist("ACGTAC", cov6)
  nb <- neighborhood("ACGTAC", 2)
  omm <- unique(unlist(lapply(nb$mm, brute_overlaps, max_extend = 3)))
  excluded <- unique(canon(c("ACGTAC", nb$mm, nb$o, omm)))
  expect_equal(sort(names(nd)), sort(setdiff(colnames(cov6), excluded)))
  expect_equal(nd, cov6["ACGTAC", names(nd)])
  expect_false("ACGTAC" %in% names(nd))
  expect_true(length(nd) > 100)
})

example_group <- data.frame(
  kmer = c("ACG", "CCG", "GCG", "TCG", "AAG", "AGG", "ATG", "ACA", "ACC",
           "ACT"),
  mismatch = c(NA, 1, 1, 1, 2, 2, 2, 3, 3, 3),
  shift = 0,
  covariability = c(1, 0.5, 0.3, -0.2, 0, 0.6, 0, 0.1, 0.1, 0.1),
  pval.adj = c(0, 0.001, 0.5, 0.001, 1, 0.001, 1, 0.5, 0.5, 0.5),
  stringsAsFactors = FALSE)

example_core <- cbind(c(0.8, 0.2, 0, 0),
                      c(0, 1 / 1.36, 0.36 / 1.36, 0),
                      c(0, 0, 1, 0))

test_that("kmer_group_to_pwm weights mismatches by squared covariability", {
  res <- chromVAR:::kmer_group_to_pwm(example_group)
  expect_equal(dim(res), c(4, 3))
  expect_equal(rownames(res), nucs)
  expect_equal(unname(res), example_core)
  expect_equal(unname(colSums(res)), rep(1, 3))
})

test_that("kmer_group_to_pwm drops mismatches above the p value threshold", {
  res <- chromVAR:::kmer_group_to_pwm(example_group, p = 0.0005)
  expect_equal(unname(res)[, 1], c(1, 0, 0, 0))
  expect_equal(unname(res)[, 2], c(0, 1, 0, 0))
})

flank_group <- function(left_cov, right_cov) {
  rbind(example_group,
        data.frame(kmer = c("AAC", "CAC", "GAC", "TAC", "CGA", "CGC", "CGG",
                            "CGT"),
                   mismatch = NA,
                   shift = c(rep(-1, 4), rep(1, 4)),
                   covariability = c(left_cov, right_cov),
                   pval.adj = ifelse(c(left_cov, right_cov) > 0, 0.001, 1),
                   stringsAsFactors = FALSE))
}

test_that("kmer_group_to_pwm adds informative flanking columns", {
  res <- chromVAR:::kmer_group_to_pwm(flank_group(c(0.8, 0.3, 0, 0),
                                                  c(0, 0, 0, 0.9)))
  expect_equal(dim(res), c(4, 5))
  left <- c(0.05 + 0.8 * 0.64 / 0.73, 0.05 + 0.8 * 0.09 / 0.73, 0.05, 0.05)
  right <- c(0.025, 0.025, 0.025, 0.925)
  expect_equal(unname(res), cbind(left, example_core, right),
               check.attributes = FALSE)
  expect_equal(unname(colSums(res)), rep(1, 5))
})

test_that("kmer_group_to_pwm trims uninformative flanking columns", {
  none <- c(0, 0, 0, 0)
  res <- chromVAR:::kmer_group_to_pwm(flank_group(none, none))
  expect_equal(unname(res), example_core)
  res <- chromVAR:::kmer_group_to_pwm(flank_group(none, c(0, 0, 0, 0.9)))
  expect_equal(dim(res), c(4, 4))
  expect_equal(unname(res)[, 4], c(0.025, 0.025, 0.025, 0.925))
  res <- chromVAR:::kmer_group_to_pwm(flank_group(c(0.8, 0.3, 0, 0), none))
  expect_equal(dim(res), c(4, 4))
  expect_equal(unname(res)[, 1],
               c(0.05 + 0.8 * 0.64 / 0.73, 0.05 + 0.8 * 0.09 / 0.73, 0.05, 0.05))
})

test_that("kmer_group_to_pwm threshold controls trimming", {
  none <- c(0, 0, 0, 0)
  res <- chromVAR:::kmer_group_to_pwm(flank_group(none, none), threshold = -1)
  expect_equal(dim(res), c(4, 5))
  expect_equal(unname(res)[, 1], rep(0.25, 4))
})

test_that("kmer_group_to_pwm ignores flanking nucleotides whose adjusted p value exceeds p", {
	g <- flank_group(c(0.8, 0.3, 0, 0), c(0, 0, 0, 0.9))
	g$pval.adj[g$shift != 0] <- 0.5
	expect_equal(unname(chromVAR:::kmer_group_to_pwm(g)), example_core)
	expect_equal(dim(chromVAR:::kmer_group_to_pwm(g, p = 0.6)), c(4, 5))
})

planted_kmer_data <- function() {
  kmers <- canonical_kmers(6)
  n <- length(kmers)
  ns <- 40
  idx <- seq_len(n * ns)
  z <- matrix(sin(idx * 1.3 + (idx %% 17) * 0.7) + cos(idx * 0.37), n, ns,
              dimnames = list(kmers, NULL))
  z <- z / sd(z)
  sig1 <- sin(seq_len(ns) * 0.9) * 3
  sig2 <- cos(seq_len(ns) * 1.7) * 3
  z["ACGTAC", ] <- z["ACGTAC", ] + sig1 * 2
  for (m in c("CCGTAC", "ACGTAA", "ACGTCC")) z[m, ] <- z[m, ] + sig1 * 1.5
  z[canon("TTGACC"), ] <- z[canon("TTGACC"), ] + sig2 * 1.8
  for (m in c("TTGACA", "TTGATC")) {
    z[canon(m), ] <- z[canon(m), ] + sig2 * 1.5
  }
  z
}

consensus <- function(motif) {
  mat <- as.matrix(motif)
  paste(rownames(mat)[apply(mat, 2, which.max)], collapse = "")
}

test_that("assembleKmers assembles a motif around a variable seed kmer", {
  z <- planted_kmer_data()
  res <- assembleKmers(SummarizedExperiment(assays = list(z = z)),
                       threshold = 2, progress = FALSE)
  expect_true(is(res, "PWMatrixList"))
  expect_equal(names(res)[1], "denovo_1")
  expect_equal(TFBSTools::name(res[[1]]), "denovo_1")
  expect_equal(TFBSTools::ID(res[[1]]), "denovo_1")
  expect_equal(TFBSTools::tags(res[[1]])$seed, "ACGTAC")
})

test_that("assembleKmers returns one motif per independent signal", {
  z <- planted_kmer_data()
  res <- assembleKmers(SummarizedExperiment(assays = list(z = z)),
                       threshold = 2, progress = FALSE)
  expect_equal(names(res), c("denovo_1", "denovo_2"))
  expect_equal(vapply(res, function(m) TFBSTools::tags(m)$seed, ""),
               c(denovo_1 = "ACGTAC", denovo_2 = canon("TTGACC")),
               ignore_attr = TRUE)
  expect_true(grepl("ACGTAC", consensus(res[[1]]), fixed = TRUE))
  expect_true(grepl(canon("TTGACC"), consensus(res[[2]]), fixed = TRUE))
  for (m in as.list(res)) {
    mat <- as.matrix(m)
    expect_equal(rownames(mat), nucs)
    expect_equal(unname(colSums(mat)), rep(1, ncol(mat)))
    expect_true(all(mat >= 0))
  }
})

test_that("assembleKmers seeds the first motif with the most variable kmer", {
  z <- planted_kmer_data()
  vars <- apply(z, 1, sd)
  res <- assembleKmers(SummarizedExperiment(assays = list(z = z)),
                       threshold = 2, progress = FALSE)
  expect_equal(TFBSTools::tags(res[[1]])$seed, names(which.max(vars)))
})

test_that("assembleKmers returns an empty list when nothing is variable", {
  z <- planted_kmer_data()
  res <- assembleKmers(SummarizedExperiment(assays = list(z = z)),
                       threshold = 100, progress = FALSE)
  expect_true(is(res, "PWMatrixList"))
  expect_equal(length(res), 0)
})

test_that("assembleKmers reports progress only when asked", {
  z <- planted_kmer_data()
  se <- SummarizedExperiment(assays = list(z = z))
  expect_output(assembleKmers(se, threshold = 2, progress = TRUE), "100%")
  expect_silent(assembleKmers(se, threshold = 2, progress = FALSE))
})

test_that("plotKmerMismatch plots one point per mismatch and seed nucleotide", {
  p <- suppressWarnings(plotKmerMismatch("ACGTAC", cov6))
  expect_s3_class(p, "ggplot")
  d <- p$data
  expect_equal(nrow(d), 18 + 6)
  expect_equal(colnames(d), c("val", "pos", "Nucleotide"))
  mm <- chromVAR:::get_mm_kmers("ACGTAC")
  cov_vals <- cov6["ACGTAC", canon(mm)]
  expect_equal(d$val[1:18], unname(ifelse(cov_vals > 0, cov_vals^2, 0)))
  expect_equal(d$pos[1:18], rep(1:6, each = 3))
  expect_equal(d$Nucleotide[1:18],
               substr(mm, rep(1:6, each = 3), rep(1:6, each = 3)))
  expect_equal(d$val[19:24], rep(1, 6))
  expect_equal(d$pos[19:24], 1:6)
  expect_equal(d$Nucleotide[19:24], strsplit("ACGTAC", "")[[1]])
  expect_true(all(d$val >= 0 & d$val <= 1))
})

test_that("plotKmerMismatch has labelled axes and one point layer", {
  p <- suppressWarnings(plotKmerMismatch("ACGTAC", cov6))
  expect_equal(p$labels$x, "Position")
  expect_equal(p$labels$y, "Shared variability with\nseed nucleotide")
  expect_equal(length(p$layers), 1)
  expect_s3_class(p$layers[[1]]$geom, "GeomPoint")
})
