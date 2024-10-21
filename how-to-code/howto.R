# excersise 1.10
A <- function(x, y) {
  print(paste0("x=", x, " ", "y=", y))
  if (y == 0) {
    return(0)
  } else if (x == 0) {
    return(2 * y)
  } else if (y == 1) {
    return(2)
  } else {
    return(A(x - 1, A(x, y - 1))) # return(A(0, 2)) # nolint
  } # call tej fun aż nie wyjdzie wynik (1,1), i jeśli y=1 to return(2), nowym y będzie 2, # nolint
  # stąd A(0,2), i jeśli x=0 to y=2y i dlatego rośnie exp.
  # dlaczego się zatrzymuje na 1024? jak działa counter?
}
A(1, 10)
A(2, 4)
A(3, 3)

# exponentiation (linear iteration)
expt <- function(b, n) {
  return(expt_iter(b, n, 1)) # ta fun jest potrzebna tylko po to, żeby ustawić product = 1 (?) # nolint
}

expt_iter <- function(b, n, product) {
  if (n == 0) {
    return(product)
  } else {
    print(product)
    return(expt_iter(b, n - 1, b * product))
  } # 2,  2-1,  2*1
  # 2,  1-1,  2*2 -> conter = 0, return(product)
}
expt(2, 2)


# 1.3.1 Procedures as Arguments
summing_fun <- function(term, a, nex_t, b) {
  if (a > b) {
    return(0)
  } else {
    return(term(a) + summing_fun(term, nex_t(a), nex_t, b))
    # uruchamia samą siebie, tylko nie z początkowym a, tylko z następnym nex_t(a) # nolint
  } # cube(1) =1               #cube, inc(1), inc, 10 # nolint
  # cube(2) =8               #cube, inc(2), inc, 10 # nolint
}

sum_cubes <- function(a, b) {
  return(summing_fun(cube, a, inc, b))
}

inc <- function(n) {
  return(n + 1)
}

cube <- function(x) {
  print(x * x * x)
  return(x * x * x)
}

sum_cubes(1, 10)
#         lower and upper bonds #nolint

identitty <- function(x) {
  return(x)
}
# jako że wykonujemy tu dowanie, które jest core funkją
# summing_fun, zastępujemy funkcję wykonującą jakieś działanie,
# funkcją identitty, która zwraca to samo, co przyjęła.
sum_integers <- function(a, b) {
  return(summing_fun(identitty, a, inc, b))
}

sum_integers(1, 10)
