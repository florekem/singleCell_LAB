library(purrr)


prof <- function(numbers) {
  nums <- c()
  alas <- map(numbers, ~ .x %% c(1:(.x - 1)))
  imap(alas[[1]], ~ ifelse(.x == 0, nums <<- c(nums, .y), FALSE))
  return(sum(nums))
}

prof(284)
prof(220)

prof(c(220, 284))

test <- function(v_of_values) {
  result <- list()
  for (i in seq_along(v_of_values)) {
    sum_of_modals <- prof(v_of_values[[i]])
    result[[i]] <- sum_of_modals
    names(result[[i]]) <- v_of_values[[i]]
  }
  return(result)
}
tutka <- test(c(220, 284))
tutka <- test(c(1:1000))
tutka
# znajdź wszystkie dzielniki liczby x
# jeśli x %% y == 0 to zapisz y.
# x <- c(1:1000) \nolint
# każdy x podziel przez y
