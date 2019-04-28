context("regioned_dna-Class")
library("SynMut")
filepath <-
    system.file("extdata", "example.fasta", package = "SynMut")
filepath.csv <-
    system.file("extdata", "target_regions.csv", package = "SynMut")
region <- read.csv(filepath.csv)
rgd.seq <- input_seq(filepath, region)

test_that("regioned_dna", {
    expect_true(class(rgd.seq)[1] == "regioned_dna")
    expect_error(input_seq(filepath, region = 1))
})

test_that("get_cu", {
    expect_true(is.matrix(get_cu(rgd.seq)))
})

test_that("get_du", {
    expect_true(is.matrix(get_du(rgd.seq)))
})

test_that("get_region", {
    expect_true(is.list(get_region(rgd.seq)))
})

test_that("get_dna", {
    expect_true(class(get_dna(rgd.seq))[1] == "DNAStringSet")
})

test_that("get_freq", {
    expect_true(is.matrix(get_freq(rgd.seq)))
})

test_that("get_rscu", {
    expect_true(is.matrix(get_rscu(rgd.seq)))
})
