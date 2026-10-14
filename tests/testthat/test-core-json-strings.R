test_that("canonical string quoting matches the established encoder byte for byte", {
 ascii <- intToUtf8(32:126,multiple=TRUE)
 pairs <- as.vector(outer(ascii,ascii,paste0))
 x <- c(ascii,pairs,intToUtf8(1:31,multiple=TRUE),NA_character_,"","</script>","<!--",
   "a/b","a\\b","\"quoted\"","a\nb","\u00e4","\u03bb","\u2028","\u2029","\U0001F4CA",intToUtf8(127))
 expected <- vapply(x,function(s) if(is.na(s)) "null" else as.character(jsonlite::toJSON(s,auto_unbox=TRUE)),"",USE.NAMES=FALSE)
 expect_identical(.ig_json_atomic_elements(x,"quote-test"),expected)
 expect_identical(.ig_json_atomic_elements(character(),"quote-test"),character())
 expect_identical(.ig_canonical_json(list(html="</script>",plain="status",empty="")),
   '{"empty":"","html":"<\\/script>","plain":"status"}')
})
