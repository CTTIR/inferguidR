th_condition_code <- function(expr) tryCatch({force(expr); NA_character_}, inferguidR_error = function(e) e$status$code)
