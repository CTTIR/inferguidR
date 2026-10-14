test_that("paragraph pagination preserves text and related values", {
  header <- c(paste("header", 1:12), "")
  analysis <- c(paste("analysis", 1:22), "Group size n1: 2", "Group size n2: 2", "")
  resources <- c("TEACHING RESOURCES", paste("reference", 1:16), "")
  limits <- c("LIMITATIONS", paste("limit", 1:7))
  lines <- c(header, analysis, resources, limits)
  pages <- .ig_report_page_chunks(lines, 50L, paragraphs = TRUE)
  expect_identical(unlist(pages, use.names = FALSE), lines)
  expect_true(all(lengths(pages) <= 50L))
  expect_length(pages, 2)
  expect_true(all(c("Group size n1: 2", "Group size n2: 2") %in% pages[[1]]))
  expect_true(all(resources %in% pages[[2]]))
})

test_that("short tails move whole paragraphs without creating a short preceding page", {
  blocks <- lapply(c(12L, 22L, 4L, 12L, 8L), function(n) c(rep("content", n - 1L), ""))
  lines <- unlist(blocks, use.names = FALSE)
  pages <- .ig_report_page_chunks(lines, 50L, paragraphs = TRUE)
  expect_identical(unlist(pages, use.names = FALSE), lines)
  expect_length(pages, 2)
  expect_gte(length(pages[[1]]), 25L)
  expect_gte(length(pages[[2]]), 20L)
})

test_that("oversized paragraphs and blank lines cannot lose or reorder content", {
  for (capacity in c(5L, 50L, 54L)) for (sizes in list(c(1, 2, 3), c(1, 101, 3), c(0, 0, 0), c(51, 2, 108), c(50, 50, 1))) {
    lines <- unlist(lapply(seq_along(sizes), function(i) c(if (sizes[i]) paste0(i, ":", seq_len(sizes[i])) else character(), "")), use.names = FALSE)
    pages <- .ig_report_page_chunks(lines, capacity, paragraphs = TRUE)
    expect_identical(unlist(pages, use.names = FALSE), lines)
    expect_true(all(lengths(pages) <= capacity))
    expect_true(all(lengths(pages) > 0L))
  }
  expect_identical(.ig_report_page_chunks(character(), 50L, paragraphs = TRUE), list())
  expect_identical(.ig_report_page_chunks("only", 50L, paragraphs = TRUE), list("only"))
})
