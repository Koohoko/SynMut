context("Mutations")
library("SynMut")
filepath <- system.file("extdata", "example.fasta", package = "SynMut")
filepath.csv <- system.file("extdata", "target_regions.csv", package = "SynMut")
region <- read.csv(filepath.csv)
rgd.seq <- input_seq(filepath, region)
rgd.seq.no.region <- input_seq(filepath)

filepath <-
    system.file("extdata", "example2.fasta", package = "SynMut")
rgd.seq.single <- input_seq(filepath)

test_that("codon_random", {
    expect_silent(seq_random(n = 1, m = 99))
    expect_silent(seq_random(n = 10, m = 30))
    expect_silent(seq_random(n = 10, m = 1:10))
    expect_silent(seq_random(
        n = 10,
        m = 100,
        no.stop.codon = TRUE
    ))
})

test_that("codon_random", {
    tmp <- codon_random(rgd.seq, n = 0.5)
    expect_true(class(tmp)[1] == "regioned_dna")
    expect_silent(codon_random(rgd.seq.single, n = 0.5))
    expect_silent(codon_random(rgd.seq.no.region, n = 0.5))

    # test numcode
    expect_silent(codon_random(rgd.seq.single, n = 0.5, numcode=2))
    ori <- get_dna(rgd.seq.single)
    mut_1 <- get_dna(codon_random(rgd.seq.single, n = 1, keep = TRUE, numcode=1)) 
    mut_2 <- get_dna(codon_random(rgd.seq.single, n = 1, keep = TRUE)) 
    mut_3 <- get_dna(codon_random(rgd.seq.single, n = 1, keep = TRUE, numcode=5)) 
    expect_true(Biostrings::translate(mut_1) == Biostrings::translate(mut_2))
    expect_false(Biostrings::translate(ori) == Biostrings::translate(mut_3))
    expect_true(Biostrings::translate(ori, genetic.code=getGeneticCode("5")) == Biostrings::translate(mut_3, genetic.code=getGeneticCode("5")))
    expect_true(all(get_cu(mut_1) == get_cu(mut_3)))
})

test_that("codon_to", {
    expect_error(codon_to(rgd.seq, max.codon = "aac", min.codon = "aac"))
    expect_error(codon_to(rgd.seq))
    expect_error(codon_to(rgd.seq, min.codon = "aaaa"))
    tmp <-
        get_cu(codon_to(rgd.seq, max.codon = "AAC")) - get_cu(rgd.seq)
    expect_true(class(tmp)[1] == "matrix")
    expect_silent(codon_to(rgd.seq.single, min.codon = "aaa"))
    expect_silent(codon_to(rgd.seq.no.region, min.codon = "aaa"))
})

test_that("dinu_to", {
    tmp <- dinu_to(rgd.seq, max.dinu = "cg")
    expect_true(class(tmp)[1] == "regioned_dna")
    tmp <- get_du(dinu_to(rgd.seq, min.dinu = "aa")) - get_du(rgd.seq)
    expect_true(class(tmp)[1] == "matrix")
    expect_silent(dinu_to(rgd.seq.single, max.dinu = "cg"))
    expect_silent(dinu_to(rgd.seq.no.region, max.dinu = "cg"))

    expect_silent(get_du(dinu_to(
        rgd.seq, max.dinu = "cg", keep = TRUE
    )))
    expect_silent(get_du(dinu_to(
        rgd.seq, min.dinu = "tg", keep = TRUE
    )))
    expect_silent(get_cu(dinu_to(
        rgd.seq.single, max.dinu = "cg", keep = TRUE
    )))
    expect_silent(get_du(dinu_to(
        rgd.seq.single, min.dinu = "tg", keep = TRUE
    )))
    expect_silent(get_cu(dinu_to(
        rgd.seq.no.region, max.dinu = "cg", keep = TRUE
    )))
    expect_silent(get_du(dinu_to(
        rgd.seq.no.region, min.dinu = "tg", keep = TRUE
    )))

    # test numcode
    ori <- get_dna(rgd.seq.single)
    mut_1 <- get_dna(dinu_to(rgd.seq.single, max.dinu = "cg", keep = TRUE, numcode=1)) 
    mut_2 <- get_dna(dinu_to(rgd.seq.single, max.dinu = "cg", keep = TRUE)) 
    mut_3 <- get_dna(dinu_to(rgd.seq.single, max.dinu = "cg", keep = TRUE, numcode=5)) 

    expect_true(Biostrings::translate(mut_1) == Biostrings::translate(mut_2))
    expect_true(Biostrings::translate(ori) == Biostrings::translate(mut_2))
    expect_true(Biostrings::translate(ori, genetic.code=getGeneticCode("5")) == Biostrings::translate(mut_3, genetic.code=getGeneticCode("5")))
    expect_true(all(get_cu(mut_1) == get_cu(mut_3)))
    expect_true(all(get_cu(ori) == get_cu(mut_3)))

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

    expect_silent(codon_mimic(rgd.seq.single, alt = target))
    expect_silent(codon_mimic(rgd.seq.no.region, alt = target))

     # test numcode
    ori <- get_dna(rgd.seq.single)
    target <- get_cu(rgd.seq)[2,]
    mut_1 <- get_dna(codon_mimic(rgd.seq.single, alt = target))
    mut_2 <- get_dna(codon_mimic(rgd.seq.single, alt = target, numcode=1))
    mut_3 <- get_dna(codon_mimic(rgd.seq.single, alt = target, numcode=5))

    idx <- which(strsplit(consensusString(c(Biostrings::translate(ori), Biostrings::translate(mut_1))), "")[[1]]=="?")

    expect_true(Biostrings::translate(ori) == Biostrings::translate(mut_1))
    expect_true(Biostrings::translate(mut_1) == Biostrings::translate(mut_2))

    # expect_false(Biostrings::translate(ori, genetic.code=getGeneticCode("1")) == Biostrings::translate(mut_3, genetic.code=getGeneticCode("1")))
    expect_true(Biostrings::translate(ori, genetic.code=getGeneticCode("5")) == Biostrings::translate(mut_3, genetic.code=getGeneticCode("5")))

})

test_that("dinu_dist", {
    expect_silent(codon_dist(codon_random(rgd.seq), rgd.seq))
    expect_silent(dinu_dist(codon_random(rgd.seq), rgd.seq))
})
