context("Mutations")
library("SynMut")
filepath <- system.file("extdata", "example.fasta", package = "SynMut")
filepath.csv <- system.file("extdata", "target_regions.csv", package = "SynMut")
region <- read.csv(filepath.csv)
rgd.seq <- input_seq(filepath, region)

filepath <- system.file("extdata", "example2.fasta", package = "SynMut")
rgd.seq2 <- input_seq(filepath)

test_that("codon_random", {
  tmp <- codon_random(rgd.seq, n = 0.5)
  expect_true(class(tmp)[1] == "regioned_dna")
  expect_silent(codon_random(rgd.seq2, n = 0.5))
})

test_that("codon_to", {
  expect_error(codon_to(rgd.seq, max.codon = "aac", min.codon = "aac"))
  expect_error(codon_to(rgd.seq))
  expect_error(codon_to(rgd.seq, min.codon = "aaaa"))
  tmp <- get_cu(codon_to(rgd.seq, max.codon = "AAC")) - get_cu(rgd.seq)
  expect_true(class(tmp)[1] == "matrix")
  expect_silent(codon_to(rgd.seq2, min.codon = "aaa"))
})

test_that("dinu_to", {
  tmp <- dinu_to(rgd.seq, max = "cg")
  expect_true(class(tmp)[1] == "regioned_dna")
  tmp <- get_du(dinu_to(rgd.seq, min = "aa")) - get_du(rgd.seq)
  expect_true(class(tmp)[1] == "matrix")
  expect_silent(dinu_to(rgd.seq2, max = "cg"))
})

test_that("codon_mimic", {
  target <- get_cu(rgd.seq)[2,]
  new <- codon_mimic(rgd.seq, alt = target)
  tmp <- get_freq(new) - get_freq(rgd.seq)
  expect_true(class(tmp)[1] == "matrix")

  target <- Biostrings::DNAStringSet("TTGAAAA-CTC-N--AAG")
  new <- codon_mimic(rgd.seq, alt = target)
  tmp <- get_cu(new) - get_cu(rgd.seq)
  expect_true(class(tmp)[1] == "matrix")

  expect_silent(codon_mimic(rgd.seq2, alt = target))
})