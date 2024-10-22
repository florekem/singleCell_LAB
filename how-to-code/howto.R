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
  if (a > b) { # ^fun     ^fun
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

# block structure, that will be removed using lambda later
pi_sum <- function(a, b) {
  pi_term <- function(x) {
    return(1.0 / (x * (x + 2)))
  }
  pi_next <- function(x) {
    return(x + 4)
  }
  return(summing_fun(pi_term, a, pi_next, b))
}
8 * pi_sum(1, 1000)

# anonymous funcions (lambda)
function(x) {
  x + 4
}

function(x) x + 4

pi_sum <- function(a, b) {
  summing_fun(
    function(x) 1.0 / (x * (x + 2)), # term ('a' is a variable for x)
    a,
    function(x) x + 4, # nex_t ('a' is a variable for x)
    b
  )
}
8 * pi_sum(1, 1000)

integral <- function(f, a, b, dx) {
  summing_fun(
    f,
    a + (dx / 2),
    function(x) x + dx,
    b
  ) * dx
}
integral(cube, 0, 1, 0.01)

(function(x, y, z) x + y + z * z)(1, 2, 3)

# f (x,y) = xa2 + yb + ab:
f <- function(x, y) {
  f_helper <- function(a, b) {
    x * (a^2) + y * b + a * b
  }
  f_helper(1 + (x * y), 1 - y) # oblicza a i b na podstawie x i y
}
f(1, 2)

f <- function(x, y) {
  # trzeba podać wartości a i b, żeby zmienić na lambdę
  a <- 1 + (x * y)
  b <- 1 - y
  (function(a, b) x * (a^2) + y * b + a * b)(a, b) # dlaczego (a,b)?
  # bo musi poznać argumenty, function(a, b) to tylko jakieś a i b,
  # a następnie treba podać jakie a i b to ma być.
}
f(1, 2)
# i jeszcze prościej:
f <- function(x, y) {
  a <- 1 + (x * y)
  b <- 1 - y
  x * (a^2) + y * b + a * b
}
f(1, 2)

# Procedures as General Methods
average <- function(a, b) {
  return((a + b) / 2.0)
}

close_enough <- function(a, b) {
  epsilon <- 0.00001
  return((b - a) < epsilon)
}

search <- function(f, neg_point, pos_point) {
  midpoint <- average(neg_point, pos_point)
  if (close_enough(neg_point, pos_point)) {
    return(midpoint)
  }

  test_value <- f(midpoint)
  if (test_value > 0) {
    return(search(f, neg_point, midpoint))
  } else if (test_value < 0) {
    return(search(f, midpoint, pos_point))
  } else {
    return(midpoint)
  }
}

# Example usage
f <- function(x) {
  return(x * x - 2)
}
result <- search(f, 0, 2)
print(result)
