test_binarize <- function() {
  test_that("binarize function works correctly", {
    data_test <- data.frame(x = c(1, 2, 3, 4, 5),
                            y = c(2, 4, 6, 8, 10))

    result <- binarize(data_test, cut_off = 3)

    expect_equal(colnames(result), c("V1", "V2", "v_bin1", "v_bin2"))
    expect_equal(result$v_bin1, c(0, 0, 0, 1, 1))
    expect_equal(result$v_bin2, c(0, 1, 1, 1, 1))
  })
}

test_pheno_frequency <- function() {
  test_that("pheno_frequency function works correctly", {
    data_test <- data.frame(v_bin1 = c(1, 1, 0, 0, 0),
                            v_bin2 = c(0, 0, 1, 1, 0))

    result <- pheno_frequency(data_test)

    expect_equal(nrow(result), 3)
    expect_equal(sum(result$freq), 5)
    expect_equal(max(result$total_bin), 1)
  })
}

test_plot_pheno <- function() {
  test_that("plot_pheno function returns a ggplot", {
    data_test <- data.frame(freq = c(2, 2, 2, 2, 1))

    result <- plot_pheno(data_test)

    expect_s3_class(result, "ggplot")
  })
}

test_describe_pheno <- function() {
  test_that("describe_pheno function works correctly", {
    data_test <- data.frame(V1 = c(1, 1, 0, 0, 0),
                            V2 = c(0, 0, 1, 1, 0),
                            freq = c(2, 2, 2, 2, 1))

    result <- describe_pheno(data_test, frequency = "freq")

    expect_equal(dim(result), c(3, 1))
    expect_equal(result[1,], 5)   # Unique Phenotypes
    expect_equal(result[2,], 2)   # N most frequent
    expect_equal(result[3,], 2)   # Median frequency
  })
}


test_common_pheno <- function() {
  test_that("common_pheno function works correctly", {
    data_test <- data.frame(V1 = c(1, 1, 0, 0, 0),
                            V2 = c(0, 0, 1, 1, 0),
                            freq = c(2, 2, 2, 2, 1))

    result <- common_pheno(data_test, n_phenotypes = 3)

    expect_equal(nrow(result), 3)  # Expecting three phenotypes
    expect_true(2 %in% result$freq) # The number 2 is the frequency of the most common phenotype
  })
}



# Run the tests
test_binarize()
test_pheno_frequency()
test_plot_pheno()
test_describe_pheno()
test_common_pheno()
