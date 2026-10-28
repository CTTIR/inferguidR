# Distribution calculator: tokenizer and parser ------------------------------------------
#
# Ownership: phase 10 utilities (group G12). Method notes:
# specs/method-notes/utilities.md, section "Distribution calculator".
#
# The calculator reads scripts in the syntax of the 'G*Power' 3.1.9.7
# calculator (manual ch. 2; conventions captured live under Wine, see the
# method notes). Input text is never evaluated as R code: it is split into
# tokens by an allowlisted lexer and parsed by a recursive-descent parser into
# a small syntax tree whose only operations are numbers, variables, the five
# arithmetic operators and calls of the functions listed in
# calculator_functions(). Evaluation is in R/calc-eval.R.
#
# Grammar (precedence from loosest to tightest; captured live):
#
#   script     := statement (SEP statement)*        SEP is ';' or a line break
#   statement  := <empty> | NAME '=' expr | expr
#   expr       := term (('+' | '-') term)*          left-associative
#   term       := power (('*' | '/') power)*        left-associative
#   power      := unary ('^' unary)*                LEFT-associative: 2^3^2 = 64
#   unary      := '-' unary | primary                binds tighter than '^': -2^2 = 4
#                                                   (no unary '+', as in G*Power)
#   primary    := NUMBER | NAME '(' expr (',' expr)* ')' | NAME | '(' expr ')'
#
# '#' starts a comment that runs to the end of the line. Names are letters
# followed by letters or digits and are case-sensitive (zCDF is not zcdf; a
# is not A). Function names cannot be used as variables. There are no
# predefined constants (pi and e are unknown variables). The argument
# separator is ','.

.ig_calc_limits <- list(max_chars = 20000L, max_depth = 200L, max_statements = 2000L)

# Tokens ----------------------------------------------------------------------------------

# Split a script into tokens. Returns list(ok = TRUE, tokens) or
# list(ok = FALSE, kind = "syntax", message, position). Every token is a list
# with `type` (num, name, op, lparen, rparen, comma, assign, sep, end), `text`,
# `value` (numbers) and `pos` (1-based character position).
.ig_calc_tokenize <- function(text) {
  chars <- strsplit(text, "", fixed = TRUE)[[1]]
  n <- length(chars)
  toks <- vector("list", 0L)
  add <- function(type, txt, pos, value = NA_real_) {
    toks[[length(toks) + 1L]] <<- list(type = type, text = txt, value = value, pos = pos)
  }
  is_digit <- function(ch) ch %in% c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9")
  # Names are letters followed by letters or digits ('_' is not part of a
  # name in G*Power: 'my_var' is read as the unknown variable 'my').
  is_alpha <- function(ch) grepl("^[A-Za-z]$", ch)
  i <- 1L
  while (i <= n) {
    ch <- chars[i]
    if (ch %in% c(" ", "\t", "\r", "\f", "\v")) {
      i <- i + 1L
    } else if (ch == "#") {
      while (i <= n && chars[i] != "\n") i <- i + 1L
    } else if (ch == "\n" || ch == ";") {
      add("sep", ch, i)
      i <- i + 1L
    } else if (is_digit(ch) || (ch == "." && i < n && is_digit(chars[i + 1L]))) {
      j <- i
      while (j <= n && is_digit(chars[j])) j <- j + 1L
      if (j <= n && chars[j] == ".") {
        j <- j + 1L
        while (j <= n && is_digit(chars[j])) j <- j + 1L
      }
      if (j <= n && chars[j] %in% c("e", "E")) {
        k <- j + 1L
        if (k <= n && chars[k] %in% c("+", "-")) k <- k + 1L
        if (k <= n && is_digit(chars[k])) {
          while (k <= n && is_digit(chars[k])) k <- k + 1L
          j <- k
        }
      }
      lit <- paste(chars[i:(j - 1L)], collapse = "")
      # The literal matches [0-9]*\.?[0-9]*([eE][+-]?[0-9]+)? by construction;
      # .ig_parse_decimal() converts it (no evaluation of input text).
      add("num", lit, i, .ig_parse_decimal(lit))
      i <- j
    } else if (is_alpha(ch)) {
      j <- i
      while (j <= n && (is_alpha(chars[j]) || is_digit(chars[j]))) j <- j + 1L
      add("name", paste(chars[i:(j - 1L)], collapse = ""), i)
      i <- j
    } else if (ch %in% c("+", "-", "*", "/", "^")) {
      add("op", ch, i)
      i <- i + 1L
    } else if (ch == "(") {
      add("lparen", ch, i)
      i <- i + 1L
    } else if (ch == ")") {
      add("rparen", ch, i)
      i <- i + 1L
    } else if (ch == ",") {
      add("comma", ch, i)
      i <- i + 1L
    } else if (ch == "=") {
      add("assign", ch, i)
      i <- i + 1L
    } else {
      return(list(ok = FALSE, kind = "syntax", position = i,
                  message = sprintf("unexpected character '%s' at position %d", ch, i)))
    }
  }
  add("end", "", n + 1L)
  list(ok = TRUE, tokens = toks)
}

# Parser ----------------------------------------------------------------------------------

# Parse a script into statements. Returns list(ok = TRUE, statements) where
# each statement is list(kind = "assign" | "expr" | "empty", name, tree,
# text, pos), or list(ok = FALSE, kind, message, position).
.ig_calc_parse <- function(text) {
  if (nchar(text) > .ig_calc_limits$max_chars) {
    return(list(ok = FALSE, kind = "syntax", position = NA_integer_,
                message = sprintf("the script is longer than %d characters", .ig_calc_limits$max_chars)))
  }
  lex <- .ig_calc_tokenize(text)
  if (!isTRUE(lex$ok)) return(lex)
  toks <- lex$tokens
  st <- new.env(parent = emptyenv())
  st$i <- 1L
  st$depth <- 0L
  peek <- function() toks[[st$i]]
  advance <- function() {
    t <- toks[[st$i]]
    st$i <- st$i + 1L
    t
  }
  fail <- function(kind, message, tok) {
    stop(structure(class = c("ig_calc_parse_error", "error", "condition"),
                   list(message = message, call = NULL, kind = kind, position = tok$pos)))
  }
  describe <- function(tok) {
    if (identical(tok$type, "end")) return("the end of the script")
    if (identical(tok$type, "sep")) return(if (tok$text == ";") "';'" else "a line break")
    sprintf("'%s'", tok$text)
  }
  enter <- function(tok) {
    st$depth <- st$depth + 1L
    if (st$depth > .ig_calc_limits$max_depth) {
      fail("syntax", sprintf("the expression is nested more than %d levels deep", .ig_calc_limits$max_depth), tok)
    }
  }
  leave <- function() st$depth <- st$depth - 1L

  parse_expr <- function() {
    enter(peek())
    on.exit(leave())
    lhs <- parse_term()
    while (identical(peek()$type, "op") && peek()$text %in% c("+", "-")) {
      op <- advance()
      rhs <- parse_term()
      lhs <- list(type = "bin", op = op$text, lhs = lhs, rhs = rhs, pos = op$pos)
    }
    lhs
  }
  parse_term <- function() {
    lhs <- parse_power()
    while (identical(peek()$type, "op") && peek()$text %in% c("*", "/")) {
      op <- advance()
      rhs <- parse_power()
      lhs <- list(type = "bin", op = op$text, lhs = lhs, rhs = rhs, pos = op$pos)
    }
    lhs
  }
  parse_power <- function() {
    lhs <- parse_unary()
    while (identical(peek()$type, "op") && peek()$text == "^") {
      op <- advance()
      rhs <- parse_unary()
      lhs <- list(type = "bin", op = "^", lhs = lhs, rhs = rhs, pos = op$pos)
    }
    lhs
  }
  parse_unary <- function() {
    tok <- peek()
    if (identical(tok$type, "op") && tok$text == "-") {
      advance()
      enter(tok)
      on.exit(leave())
      arg <- parse_unary()
      return(list(type = "neg", arg = arg, pos = tok$pos))
    }
    if (identical(tok$type, "op") && tok$text == "+") {
      # G*Power 3.1.9.7 refuses a unary plus ('+3' gives an error box; captured
      # live, SUP-G12-UT-CALC-P07); the calculator accepts the same language.
      fail("syntax", "a unary '+' is not part of the calculator syntax", tok)
    }
    parse_primary()
  }
  parse_primary <- function() {
    tok <- peek()
    if (identical(tok$type, "num")) {
      advance()
      return(list(type = "num", value = tok$value, text = tok$text, pos = tok$pos))
    }
    if (identical(tok$type, "name")) {
      advance()
      if (identical(peek()$type, "lparen")) {
        open <- advance()
        args <- list()
        if (identical(peek()$type, "rparen")) {
          fail("syntax", sprintf("missing argument for function %s", tok$text), peek())
        }
        repeat {
          if (peek()$type %in% c("comma", "rparen")) {
            fail("syntax", sprintf("missing argument for function %s", tok$text), peek())
          }
          args[[length(args) + 1L]] <- parse_expr()
          nxt <- peek()
          if (identical(nxt$type, "comma")) {
            advance()
            next
          }
          if (identical(nxt$type, "rparen")) {
            advance()
            break
          }
          if (nxt$type %in% c("end", "sep")) fail("syntax", "unbalanced parentheses", open)
          fail("syntax", sprintf("expected ',' or ')' in the arguments of %s, found %s", tok$text, describe(nxt)), nxt)
        }
        return(list(type = "call", fun = tok$text, args = args, pos = tok$pos))
      }
      return(list(type = "var", name = tok$text, pos = tok$pos))
    }
    if (identical(tok$type, "lparen")) {
      advance()
      inner <- parse_expr()
      if (!identical(peek()$type, "rparen")) {
        if (peek()$type %in% c("end", "sep")) fail("syntax", "unbalanced parentheses", tok)
        fail("syntax", sprintf("expected ')' but found %s", describe(peek())), peek())
      }
      advance()
      return(list(type = "group", arg = inner, pos = tok$pos))
    }
    fail("syntax", sprintf("missing operand at %s", describe(tok)), tok)
  }

  statements <- list()
  out <- tryCatch({
    repeat {
      start <- peek()
      if (start$type %in% c("sep", "end")) {
        statements[[length(statements) + 1L]] <- list(kind = "empty", name = NA_character_, tree = NULL,
                                                      pos = start$pos)
      } else if (identical(start$type, "name") && identical(toks[[st$i + 1L]]$type, "assign")) {
        advance()
        advance()
        if (peek()$type %in% c("sep", "end")) fail("syntax", sprintf("missing value in the assignment to %s", start$text), peek())
        tree <- parse_expr()
        statements[[length(statements) + 1L]] <- list(kind = "assign", name = start$text, tree = tree,
                                                      pos = start$pos)
      } else {
        tree <- parse_expr()
        statements[[length(statements) + 1L]] <- list(kind = "expr", name = NA_character_, tree = tree,
                                                      pos = start$pos)
      }
      nxt <- peek()
      if (identical(nxt$type, "end")) break
      if (identical(nxt$type, "sep")) {
        advance()
        next
      }
      if (identical(nxt$type, "rparen")) fail("syntax", "unbalanced parentheses", nxt)
      if (identical(nxt$type, "assign")) fail("syntax", "'=' is only allowed after a variable name at the start of a statement", nxt)
      fail("syntax", sprintf("missing operator before %s", describe(nxt)), nxt)
    }
    if (length(statements) > .ig_calc_limits$max_statements) {
      fail("syntax", sprintf("more than %d statements", .ig_calc_limits$max_statements), toks[[1]])
    }
    list(ok = TRUE, statements = statements)
  }, ig_calc_parse_error = function(e) {
    list(ok = FALSE, kind = e$kind, message = conditionMessage(e), position = e$position)
  })
  out
}

# Text of a syntax tree (canonical form, fully parenthesised where the
# grammar needs it). Used in notes and in the statement table.
.ig_calc_deparse <- function(node) {
  switch(node$type,
    num = node$text,
    var = node$name,
    group = paste0("(", .ig_calc_deparse(node$arg), ")"),
    neg = paste0("-", .ig_calc_deparse(node$arg)),
    bin = paste0(.ig_calc_deparse(node$lhs), node$op, .ig_calc_deparse(node$rhs)),
    call = paste0(node$fun, "(", paste(vapply(node$args, .ig_calc_deparse, ""), collapse = ","), ")"),
    stop("unknown node type", call. = FALSE)
  )
}
