library(sloop)
library(lobstr)

sloop::otype()
typeof()
attributes()
ftype()
sloop:::ftype
ftype(t.test)
# [1] "S3"      "generic"
ftype(t.data.frame)
# [1] "S3"     "method"
s3_get_method(t.data.frame)
s3_get_method(weighted.mean.Date)
t.test

# S3
# constructors
# • Be called new_myclass(). • Have one argument for the base object, and one for each attribute. • Check the type of the base object and the types of each attribute. #nolint
new_Date <- function(x = double()) {
  stopifnot(is.double(x))
  structure(x, class = "Date")
}
new_Date(0)

# validators
# Rather than encumbering the constructor with complicated checks, it’s better to put them in a separate function. Doing so allows you to cheaply create new objects when you know that the values are correct, and easily re-use the checks in other places. #nolint

# helpers

my_new_generic <- function(x) {
  UseMethod("my_new_generic")
}


s3_dispatch(mean(2.0))

f <- factor(c("a", "b", "c"))
f
class(f)
typeof(f)
attributes(f)
unclass(f)
unclass

#
x <- c(1, 2, 3)
obj_addr(x)
# [1] "0x55d5c3f7f358"
y <- x
obj_addr(y)
# [1] "0x55d5c3f7f358"
# name binds to the same object

y[[3]] <- 4
obj_addr(y)
# [1] "0x55d5c3fc3e18"
# new object as a result of modification

# If you modify y again, it won’t get copied.
# That’s because the new object now only has a single name
# bound to it, so R applies modify-in-place optimisation.
y[[3]] <- 5
obj_addr(y)

x <- list(1:10)

?rlang::env()

# You can determine the type of a vector with typeof()
# and its length with length().
typeof(x)
# [1] "double"
int_vect <- c(1L, 2L)
typeof(int_vect)
# [1] "integer"

# coercion
# When you attempt to combine different types they will be coerced
# in a fixed order: character → double → integer → logical.
str(c("a", 1))
#> chr [1:2] "a" "1"

?structure()

x <- c(FALSE, FALSE, TRUE)
as.numeric(x) #> [1] 0 0 1
# Total number of TRUEs
sum(x)
#> [1] 1
# Proportion that are TRUE
mean(x)
#> [1] 0.333

otype(1:10)

sloop::s3_class(matrix(1:4, nrow = 2))
sloop::otype(matrix(1:4, nrow = 2))

# Attributes
str(1:3) # 1d vector
dim(str(1:3)) # vector w/o a dim attribute has NULL dims

str(matrix(1:3)) # column vector
#> int [1:3, 1] 1 2 3

str(matrix(1:3, ncol = 1)) # column vector
#> int [1:3, 1] 1 2 3

str(matrix(1:3, nrow = 1)) # row vector
#> int [1, 1:3] 1 2 3

str(array(1:3, 3)) # "array" vector
#> int [1:3(1d)] 1 2 3

# vector could be given factor attribute,
# and converted into S3
sex_char <- c("m", "m", "m")
s3_class(sex_char)
#> [1] "character"
attributes(sex_char)
#> NULL
sex_factor <- factor(sex_char, levels = c("m", "f"))
s3_class(sex_factor)
# [1] "factor"
attributes(sex_factor)
# $levels
# [1] "m" "f"

# $class
# [1] "factor"
table(sex_char)
#> sex_char
#> m
#> 3
# when you tabulate a factor you’ll get counts
# of all categories, even unobserved ones
table(sex_factor)
#> sex_factor
#> m f
#> 3 0

# lists
# c() will combine several lists into one.
l4 <- list(list(1, 2), c(3, 4))
str(l4)
# List of 2
#  $ :List of 2
#   ..$ : num 1
#   ..$ : num 2
#  $ : num [1:2] 3 4

# If given a combination of atomic vectors and lists,
# c() will coerce the vectors to lists before combining them.
l5 <- c(list(1, 2), c(3, 4))
str(l5)
# List of 4
#  $ : num 1
#  $ : num 2
#  $ : num 3
#  $ : num 4

# With atomic vectors, the dimension attribute
# is commonly used to create matrices.
x <- c(1, 2, 3, 4)
dim(x) <- c(2, 2)
x
#      [,1] [,2]
# [1,]    1    3
# [2,]    2    4
dim(x)
# [1] 2 2
str(x)
# num
sloop::otype(x)
# "base"
attributes(x)
# $dim
# [1] 2 2
typeof(x)
# "double"

# With lists, the dimension attribute can be used
# to create list-matrices or list-arrays:
l <- list(1:3, "a", TRUE, 1.0)
l
dim(l) <- c(2, 2)
#      [,1]      [,2]
# [1,] integer,3 TRUE
# [2,] "a"       1
sloop::otype(l)
# "base"
attributes(l)
# $dim
# [1] 2 2
typeof(l)
# "list"

# A data frame is a named list of vectors with attributes
# for (column) names, row.names11, and its class, “data.frame”:
df1 <- data.frame(x = 5:7, y = letters[1:3])
df1
#   x y
# 1 5 a
# 2 6 b
# 3 7 c
attributes(df1)
# $names (colnames)
# [1] "x" "y"

# $class
# [1] "data.frame"

# $row.names
# [1] 1 2 3
typeof(df1)
# "list"
sloop::otype(df1)
# "S3"

df1 <- data.frame(
  x = 1:3,
  y = c("a", "b", "c"),
  stringsAsFactors = FALSE
)

library(tibble)
df2 <- tibble(
  x = 5:7,
  y = c("a", "b", "c")
)
str(df2)
# tibble [3 × 2] (S3: tbl_df/tbl/data.frame)
#  $ x: int [1:3] 5 6 7
#  $ y: chr [1:3] "a" "b" "c"
sloop::otype(df2)
# "S3"
df3 <- data.frame(
  age = c(35, 27, 18),
  hair = c("blond", "brown", "black"),
  row.names = c("Bob", "Susan", "Sam")
)
df3
# you can use rownames to subset rows:
df3["Bob", ]

# tibbles do not support rownames.
# to use rownames, as a new column"
df4 <- as_tibble(df3, rownames = "name") # create new column "name"
#   name    age hair
#   <chr> <dbl> <chr>
# 1 Bob      35 blond
# 2 Susan    27 brown
# 3 Sam      18 black

# subsetting of df will return a vector if subsetting one var
df3[, "hair"]
# [1] "blond" "brown" "black"
str(df3[, "hair"])
# chr [1:3]

# subsetting tibble always returns a tibble
df4[, "hair"]
# A tibble: 3 × 1
#   hair
#   <chr>
# 1 blond
# 2 brown
# 3 black

# IF YOU WANT A SINGLE COLUMN USE [[]]:
df3[["hair"]]
str(df3[["hair"]])
# chr
# however it will retrn a vectro even in tibble!
df4[["hair"]]
str(df4[["hair"]])

is_tibble(df4)

# List columns
# http://r4ds.had.co.nz/many-models.html
df <- data.frame(x = 1:3)
df$y <- list(1:2, 1:3, 1:4)
df
#   x          y
# 1 1       1, 2
# 2 2    1, 2, 3
# 3 3 1, 2, 3, 4

# to create df, I() must be used
# I() is short for identity and is often used to indicate
# that an input should be left as is,
# and not automatically transformed.
data.frame(
  x = 1:3,
  y = I(list(1:2, 1:3, 1:4))
)

# this is easier with tibbles:
# (and they will be printed tidily)
tibble(
  x = 1:3, y = list(1:2, 1:3, 1:4)
)
# A tibble: 3 × 2
#       x y
#   <int> <list>
# 1     1 <int [2]>
# 2     2 <int [3]>
# 3     3 <int [4]>

as.matrix(df3)
#       age  hair
# Bob   "35" "blond"
# Susan "27" "brown"
# Sam   "18" "black"
data.matrix(df3)
#       age hair
# Bob    35    2
# Susan  27    3
# Sam    18    1
data.matrix()
# Return the matrix obtained by converting all the
# variables in a data frame to numeric mode and
# then binding them together as the columns of a matrix.
# Factors and ordered factors are replaced by their internal codes.

a <- matrix(1:9, nrow = 3)
a
#      [,1] [,2] [,3]
# [1,]    1    4    7
# [2,]    2    5    8
# [3,]    3    6    9

df <- data.frame(x = 1:3, y = 3:1, z = letters[1:3])
df
df[df$x == 2, ]
df[[2]] == df[, 2]

df[c(1, 3), ]

# subsetting single column of df like a list:
list <- list("x" = c(1.1, 1.2), "y" = c(2.1, 2.2), z = c(3.1, 3.2))
list
list["x"]
list[1]
list[1:2]
list[1]
list[[1]]

str(df["x"]) # do not simplify
# subsettign single column of df like a matrix:
str(df[, "x"]) # simplifies
# int
str(df[["x"]]) # this also simplifies
# int

df <- tibble::tibble(x = 1:3, y = 3:1, z = letters[1:3])
# Subsetting a tibble with [ always returns a tibble:
str(df["x"]) # tibble
str(df[, "x"]) # tibble
# however [[]] does simplify
str(df[["x"]]) # int

a <- matrix(1:4, nrow = 2)
a
a[1, , drop = FALSE]
#      [,1] [,2]
# [1,]    1    3

View(mtcars)
mtcars[mtcars$cyl == 4, ]
mtcars[-c(1:4), ]

x <- outer(1:5, 1:5, FUN = "*")
x
x[upper.tri(x)]

df[1, 1] <- NA
#       x     y z
#   <int> <int> <chr>
# 1    NA     3 a
# 2     2     2 b
# 3     3     1 c
df[is.na(df)] <- 0 # ciekawe
#       x     y z
#   <int> <int> <chr>
# 1     0     3 a
# 2     2     2 b
# 3     3     1 c

x <- list(1:3, "a", 4:6)
x
# [[1]]
# [1] 1 2 3

# [[2]]
# [1] "a"

# [[3]]
# [1] 4 5 6

x[1] # extract a smaller train (with choochoo)
# [[1]] # choochoo
# [1] 1 2 3 # carriage
x[[1]] # extract only carriage
# [1] 1 2 3 # carriage
x[[1]][[2]] # recursive subsetting
# 2

# subsetting with nothing can be useful with assignment,
# because it prevents the structure of the original object:
mtcars[] <- lapply(mtcars, as.integer)
is.data.frame(mtcars)
#> [1] TRUE
mtcars <- lapply(mtcars, as.integer)
is.data.frame(mtcars)
#> [1] FALSE


# lookup tables (character subsetting)
x <- c("m", "f", "u", "f", "f", "m", "m")
lookup <- c(m = "Male", f = "Female", u = NA)
lookup
lookup[x]
#      m        f        u        f        f        m        m
# "Male" "Female"       NA "Female" "Female"   "Male"   "Male

# if you don't want names in the result:
unname(lookup[x])
# [1] "Male"   "Female" NA       "Female" "Female" "Male"   "Male"

# matching and merging by hand (integer subsetting)
grades <- c(1, 2, 2, 3, 1)
info <- data.frame(
  grade = 3:1,
  desc = c("Excellent", "Good", "Poor"),
  fail = c(F, F, T)
)
info
#   grade      desc  fail
# 1     3 Excellent FALSE
# 2     2      Good FALSE
# 3     1      Poor  TRUE

# Then, let’s say we want to duplicate the info table so that we have a row for each value in grades. #nolint
grades
# [1] 1 2 2 3 1
info$grade
# [1] 3 2 1
# An elegant way to do this is by combining match() and integer subsetting
# (match(needles, haystack) returns the position
# _where_ each needle is found in the haystack).
id <- match(grades, info$grade)
# _where_ each grade if found in the info$grade
id
# [1] 3 2 2 1 3 # pozycja 3, 2 , itd w info$grade
info[id, ]
#  grade      desc  fail
# 3       1      Poor  TRUE
# 2       2      Good FALSE
# 2.1     2      Good FALSE
# 1       3 Excellent FALSE
# 3.1     1      Poor  TRUE
info[id, ]$fail
# [1]  TRUE FALSE FALSE FALSE  TRUE

# if matching on multiple colums, its better to use merge()
# or dplyr::left_join()

# random samples and bootstraps (integer subsetting)
df <- data.frame(
  x = c(1, 2, 3, 1, 2), y = 5:1, z = letters[1:5]
)
# randomly reorder
sample(5)
# [1] 3 5 1 4 2
df[sample(nrow(df)), ]
df[c(3, 5, 1, 4, 2), ]
# select 3 random rows:
sample(5, 3)
# [1] 1 4 5
df[sample(nrow(df), 3), ]

# select 6 bootstrap replicates:
sample(4, 3, replace = TRUE)
# [1] 2 4 4
df[sample(nrow(df), 6, replace = TRUE), ]

# dokończyć !!!!!!!!!!!1


# 5. Control flow
# 5.2 Choices
# vectorized if
# handling _vectors_ [c(1:10)] of values is the job of ifelse()
# ifelse(cond, true, false)
# use ifelse() only when the 'yes' and 'no' vectors are the same type

# another vectorised equivalent is dplyr::caes_when()
x <- c(1:10)
dplyr::case_when(
  x %% 35 == 0 ~ "fizz buzz",
  x %% 5 == 0 ~ "fizz",
  x %% 7 == 0 ~ "buzz",
  is.na(x) ~ "???",
  .default = as.character(x)
)
#  [1] "1"    "2"    "3"    "4"    "fizz" "6"    "buzz" "8"    "9"    "fizz"
dplyr::case_when(
  x %% 35 == 0 ~ "fizz buzz",
  x %% 5 == 0 ~ "fizz",
  x %% 7 == 0 ~ "buzz",
  is.na(x) ~ "???",
  TRUE ~ as.character(x) # nie rozumiem tego zapisu, w help jest .default
)
# [1] "1"    "2"    "3"    "4"    "fizz" "6"    "buzz" "8"    "9"    "fizz"
# wynik jest ten sam

# Switch()
# I recommend usign switch() only with character inputs!!!
x_option <- function(x) {
  switch(x,
    a = "option 1", # nie musi być "a", jeśli daje string? dlaczego?
    b = "option 2",
    c = "option 3",
    stop("Invalid `x` value") # else
  )
}
x_option("a")
# The last component of a switch() should always throw an error, otherwise unmatched inputs will invisibly return NULL  # nolint
x_option("d")
# Error in x_option("d") : Invalid `x` value
(switch("c",
  a = 1,
  b = 2
))
# NULL
(switch("c",
  a = 1,
  b = 2,
  stop("invalid")
))
# Error: invalid

# Loops
means <- c(1, 50, 20)
out <- vector("list", length(means))
out
?vector
# vector(mode = "logical", length = 0) #nolint
# as.vector(x, mode = "any") # nolint
# is.vector(x, mode = "any") # nolint
for (i in seq_along(means)) {
  out[[i]] <- rnorm(10, means[[i]])
}
out

# when looping though S3 vectors, it is best to use [[
# as it prevents from stripping the attributes:
xs <- as.Date(c("2020-01-01", "2010-01-01"))
for (i in seq_along(xs)) {
  print(xs[[i]])
}

x <- numeric()
x
# numeric(0)
length(x)
# 0
out <- vector("list", length(x))
out
# list()
length(out)
# 0
out2 <- vector("list")
out2
# list()
length(out2)
# 0
out3 <- list()
out3
# list()

for (i in 1:length(x)) {
  out[i] <- x[i]^2
}
out
# zwraca NA
1:length(0)
# 1
for (i in seq_len(x)) {
  out[i] <- x[i]^2
}
# error

# 6. Functions

lapply(mtcars, function(x) length(unique(x)))

Filter(function(x) !is.numeric(x), mtcars)
# is not the same as filter() (lower case)

integrate(function(x) sin(x)^2, 0, pi)

# You can put functions in a list
funs <- list(
  half = function(x) x / 2,
  double = function(x) x * 2
)
funs$double(10)

# if you have the arguments already in data structure,
# you can use do.call(), with two arguments: function to call,
# and a list containing the function arguments
args <- list(1:10, na.rm = TRUE)
args
do.call(mean, args)

g11 <- function() {
  if (!exists("a")) {
    a <- 1
  } else {
    a <- a + 1
  }
  a
}
g11()

f <- function(x) {
  f <- function(x) {
    f <- function() {
      x^2
    }
    f() + 1
  }
  f(x) * 2
}
f(10)

h01 <- function(x) {
  10
  # x is not present in a function, x is not evaluated so no error, only 10 is printed # nolint
}
h01(stop("error"))
# 10

x <- 20
double <- function(x) {
  message("Calculating...")
  x * 2
}

h03 <- function(x) {
  c(x, x)
}
h03(double(x))


x_ok <- function(x) {
  !is.null(x) && length(x) == 1 && x > 0
}
a <- x_ok(1:3)
a
length(a)


f1 <- function(x = {
                 y <- 1
                 2
               }, y = 0) {
  c(x, y)
}
f1()

x <- {
  y <- 1
  2
}
x

?sum

browser()
plot(1:10, col = "red", pch = 20, xlab = "x", col.lab = "blue")

add <- function(x, y) x + y
a <- list(1:3, 4:5)
a
lapply(a, add, 3)
#          x   y   whyyyy?

add <- function(x) x + 3
lapply(a, add)

# the same will be obtained using:
lapply(list(1:3, 4:5), `+`, 3)


library(rlang)
e3 <- env(x = 1, y = 2)
e3$x
e3$z <- 3
e3[["z"]]

# Environment recurssion
fget <- function(name, env = caller_env(), inherits = TRUE) {
  if (identical(env, emptyenv())) {
    # base case
    stop(name, " is not a function", call. = FALSE)
  } else if (env_has(env, name) && is.function(env_get(env, name))) {
    # success case
    return("is a function")
  } else if (inherits == TRUE) {
    # recursive case
    fget(name, env_parent(env))
  } else {
    stop(name, " is not a function, but inherits == FALSE", call. = FALSE)
  }
}

e4 <- env(x = 1, y = 2, add2 = function(x = 1, y = 2) x + y)
e5 <- env(e4, x = 2, y = 4, add2 = 5, add3 = function(x, y) x + y)


fget("add2", e5)
fget("add2", e5, inherits = FALSE)

is.function(e4[["add2"]])

# conditions
# this is a function wrapper, that wraps it before args are send to
# the main funcion. It is also possible to wrap an output of the funtion
# and add additional logic to it
file_remove_strict <- function(path) {
  if (!file.exists(path)) {
    stop("could not find a file", call. = FALSE)
  }
  file.remove(path)
}
saveRDS(mtcars, "mtcars.rds")
file_remove_strict("mtcars.rds")

tryCatch(
  error = function(cnd) 10, # no error, keep going
  1 + 1
)
# but:
tryCatch(
  message = function(cnd) "There", # message, not error, will be displayed
  { # and this block will be never shown, because tryCatch() is an exit handler
    message("Here")
    stop("This code is never run!")
  }
)

catch_cnd <- function(expr) {
  tryCatch(
    condition = function(cnd) cnd,
    {
      force(expr)
      return(NULL)
    }
  )
}

catch_cnd(stop("go"))
catch_cnd(4 + 4)

catch_cnd(stop("An error"))
catch_cnd(abort("An error"))

f3 <- function(x) {
  tryCatch(
    error = function(cnd) NA,
    log(x)
  )
}
f3(4)


# Functionals
randomise <- function(f) f(runif(1e3))

randomise <- function(f) {
  f(runif(1000))
}

randomise(mean)
randomise(sum)

library(purrr)

# this is how purrr::map works:
simple_map <- function(x, f, ...) {
  out <- vector("list", length(x))
  for (i in seq_along(x)) {
    out[[i]] <- f(x[[i]], ...)
  }
  out
}
#            x        f
simple_map(mtcars, typeof)
typeof(mtcars[[1]])

x <- mac(1:3, function(x) runif(2))
x <- map(1:3, ~ runif(2))
str(x)
# List of 3
#  $ : num [1:2] 0.972 0.584
#  $ : num [1:2] 0.816 0.321
#  $ : num [1:2] 0.188 0.482

x <- list(1:5, c(1:10, NA))
map_dbl(x, function(x) mean(x, na.rm = TRUE))
# OR
map_dbl(x, ~ mean(.x, na.rm = TRUE))
# OR (as map() functions pass ... along)
map_dbl(x, mean, na.rm = TRUE)

x <- c(3, 0, 0, 0)
map_dbl(x, function(x, y) x + y, 1) # y = 1 as it is passed as a 2nd arg
#                                                       to the funcion.
# [1] 4 1 1 1
# for example:
map(x, mean, 0.1)
# will call:
map(x[[1]], 0.1) # and second arg for mean is trim = 0.1
# in the same way function(x, y) = function(arg1=x, arg2=y=1)

View(mtcars) # numeric df
str(mtcars)
map_dbl(mtcars, sd)
map_dbl(mtcars, ~ sd(.x))
map_dbl(mtcars, function(x) sd(x))

penguins <- palmerpenguins::penguins
num_only <- map_lgl(penguins, is.numeric)
penguins[num_only]
map_dbl(penguins[num_only], sd, na.rm = TRUE)

factors_ping <- map_lgl(penguins, is.factor)
factors_ping <- map_lgl(penguins, ~ is.factor(.x))
factors_ping
map_dbl(penguins[factors_ping], ~ length(levels(.x)))
# dlaczego raz działa samo sd, a innym razem trzeba przez function??
# bo wystarczy pojedyncza ewaluacja? Bez konieczności loopowania?

trials <- map(1:100, ~ t.test(rpois(10, 10), rpois(7, 10)))
str(trials)
# select elements by name
pvals <- map_dbl(trials, "p.value")
pvals

x <- list(
  list(1, c(3, 9)),
  list(c(3, 6), 7, c(4, 7, 6))
)
triple <- function(x) x * 3
map(x, map, ~ triple(.x))
map(x, map, triple)

str(mtcars)

# this base function splits df into 3 df based on no. of cylinders:
by_cyl <- split(mtcars, mtcars$cyl)
str(by_cyl)
# list of 3
by_cyl <- by_cyl |>
  map(~ lm(mpg ~ wt, data = .x))
str(by_cyl)
# 3 lists of linear models
by_cyl <- split(mtcars, mtcars$cyl)
by_cyl <- by_cyl |>
  map(~ lm(mpg ~ wt, data = .x)) |>
  map("coefficients")
by_cyl
# (Intercept)          wt
#   39.571196   -5.647025
by_cyl <- split(mtcars, mtcars$cyl)
by_cyl <- by_cyl |>
  map(~ lm(mpg ~ wt, data = .x)) |> # calculate lm
  map("coefficients") |> # extract value
  map_dbl(2) # extract wt (only 2nd element of lists) select by possition
by_cyl
#  4         6         8
# -5.647025 -2.780106 -2.192438

x <- list(c(12, 22))
x
# [1] 12 22
map_dbl(x, 2) # select by possition
# 22

trims <- c(0, 0.1, 0.2, 0.5)
x <- rcauchy(1000)
x
map_dbl(trims, ~ mean(x, trim = .x))
# in this case x in constant (we want to check different ammounts of
# trimming on a constant mean).
# this is confusing, try to make it a little bit clearer:
map_dbl(trims, function(trim) mean(x, trim = trim))
# or using pmap
pmap_dbl(list(trim = trims), mean, x = x)

modify(mtcars, 3)
?mtcars

temp <- tempfile()
dir.create(temp)
cyls <- split(mtcars, mtcars$cyl)
paths <- file.path(temp, paste0("cyl-", names(cyls), ".csv"))
walk2(cyls, paths, write.csv)
dir(temp)
# rewrite this ^^ with iwalk() instead of walk2()
temp2 <- tempdir()
dir(temp2)
iwalk(cyls, function(x, y) write.csv(x, paste0(temp, "/", y, ".csv")))
dir(temp)
# this was my interpretation of iwalk, but HWickham suggests:
# we can store path in variable names! and use that to save
# will I ever have this kind of imagination? anyway:
temp <- tempfile()
dir.create(temp)
cyls <- split(mtcars, mtcars$cyl)
names(cyls) <- file.path(temp, paste0("cyl-", names(cyls), ".csv"))
cyls
iwalk(cyls, ~ write.csv(.x, .y))
#                       df, index = names = path = genius!

funs_list <- list(
  disp = function(x) x * 0.0163871,
  am = function(x) factor(x, labels = c("auto", "manual"))
)
nm <- names(funs_list)
nm
mtcars[nm] <- map2(funs_list, mtcars[nm], function(f, var) f(var))
# which is (probably) the same as (but mayby easier to understand):
mtcars[nm] <- map2(funs_list, mtcars[nm], ~ .x(.y))
# co jest mega dziwne?
# w miejsce wektora jest lista funcji? ok, lista jest wektorem,
# i jako taka jest traktowana przez map2. wektorem jest też mtcars[mn],
# czyli się zgadza, mamy dwa wektory, które następnie przekazywane są
# do funkcji f(var), czyli z listy (wektora) stają się już funkcją,
# w których modyfikowane są kolumny disp i am. var przekazywane jest do
# x funcji disp i am?
# opis HWickhama:
# the list of the 2 functions (trans) and the 2 appropriately
# selected data frame columns (mtcars[nm]) are supplied
# to map2().
# map2() creates an anonymous function (f(var))
# which *applies the functions to the variables* when map2()
# iterates over their (similar) indices. On the left-hand side,
# the respective 2 elements of mtcars are being replaced
# by their new transformations.
rm(mtcars)

mtcars[nm]

# reduce
reduce(1:6, sum)
accumulate(1:6, sum)
# sum(1) -> 1 #nolint
# sum(1, 2) -> 3 #nolint
# sum(3, 3) -> 6 #nolint
# sum(6, 4) -> 10 #nolint
# sum(10, 5) -> 15 #nolint
# sum(15, 6) -> 21 #nolint

library(rlang)
library(ggplot2)
library(scales)

power1 <- function(exp) {
  force(exp)
  function(x) {
    x^exp
  }
}
square <- power1(2)
cube <- power1(3)

square
env_print(square)

?ecdf

chatty <- function(f) {
  force(f)
  function(x, ...) {
    res <- f(x, ...)
    cat("Processing ", x, "\n", sep = "")
    res
  }
}

f <- function(x) x^2
s <- c(3, 2, 1)
purrr::map_dbl(s, chatty(f))

urls <- c(
  "adv-r" = "https://adv-r.hadley.nz",
  "r4ds" = "http://r4ds.had.co.nz/" #
)

path <- paste(tempdir(), names(urls), ".html")

delay_by <- function(f, amount) {
  force(f)
  force(amount)

  function(...) {
    Sys.sleep(amount)
    f(...)
  }
}

dot_every <- function(f, n) {
  force(f)
  force(n)

  i <- 0
  function(...) {
    i <<- i + 1
    if (i %% n == 0) cat(".")
    f(...)
  }
}

purrr::walk2(urls, path, dot_every(delay_by(download.file, 1), 1))

purrr::walk2(urls, path, download.file |> dot_every(10) |> delay_by(10))
# dlaczego download.file nie jest wykonywany od razu i dopiero wynik jest
# pipowany dalej, tylko razem z urls i path jest pipowany do dot_every i
# delay_by?


?match.arg

library(Seurat)
Seurat:::Read10X

function(
    data.dir, gene.column = 2, cell.column = 1, unique.features = TRUE,
    strip.suffix = FALSE) {
  full.data <- list()
  has_dt <- requireNamespace("data.table", quietly = TRUE) &&
    requireNamespace("R.utils", quietly = TRUE)
  for (i in seq_along(along.with = data.dir)) {
    run <- data.dir[i]
    if (!dir.exists(paths = run)) {
      stop("Directory provided does not exist")
    }
    barcode.loc <- file.path(run, "barcodes.tsv")
    gene.loc <- file.path(run, "genes.tsv")
    features.loc <- file.path(run, "features.tsv.gz")
    matrix.loc <- file.path(run, "matrix.mtx")
    pre_ver_3 <- file.exists(gene.loc)
    if (!pre_ver_3) {
      addgz <- function(s) {
        return(paste0(s, ".gz"))
      }
      barcode.loc <- addgz(s = barcode.loc)
      matrix.loc <- addgz(s = matrix.loc)
    }
    if (!file.exists(barcode.loc)) {
      stop("Barcode file missing. Expecting ", basename(path = barcode.loc))
    }
    if (!pre_ver_3 && !file.exists(features.loc)) {
      stop(
        "Gene name or features file missing. Expecting ",
        basename(path = features.loc)
      )
    }
    if (!file.exists(matrix.loc)) {
      stop(
        "Expression matrix file missing. Expecting ",
        basename(path = matrix.loc)
      )
    }
    data <- readMM(file = matrix.loc)
    if (has_dt) {
      cell.barcodes <- as.data.frame(data.table::fread(barcode.loc,
        header = FALSE
      ))
    } else {
      cell.barcodes <- read.table(
        file = barcode.loc, header = FALSE,
        sep = "\t", row.names = NULL
      )
    }
    if (ncol(x = cell.barcodes) > 1) {
      cell.names <- cell.barcodes[, cell.column]
    } else {
      cell.names <- readLines(con = barcode.loc)
    }
    if (all(grepl(pattern = "\\-1$", x = cell.names)) & strip.suffix) {
      cell.names <- as.vector(x = as.character(x = sapply(
        X = cell.names,
        FUN = ExtractField, field = 1, delim = "-"
      )))
    }
    if (is.null(x = names(x = data.dir))) {
      if (length(x = data.dir) < 2) {
        colnames(x = data) <- cell.names
      } else {
        colnames(x = data) <- paste0(i, "_", cell.names)
      }
    } else {
      colnames(x = data) <- paste0(
        names(x = data.dir)[i],
        "_", cell.names
      )
    }
    if (has_dt) {
      feature.names <- as.data.frame(data.table::fread(ifelse(test = pre_ver_3,
        yes = gene.loc, no = features.loc
      ), header = FALSE))
    } else {
      feature.names <- read.delim(
        file = ifelse(test = pre_ver_3,
          yes = gene.loc, no = features.loc
        ), header = FALSE,
        stringsAsFactors = FALSE
      )
    }
    if (any(is.na(x = feature.names[, gene.column]))) {
      warning("Some features names are NA. Replacing NA names with ID from the opposite column requested",
        call. = FALSE, immediate. = TRUE
      )
      na.features <- which(x = is.na(x = feature.names[
        ,
        gene.column
      ]))
      replacement.column <- ifelse(test = gene.column ==
        2, yes = 1, no = 2)
      feature.names[na.features, gene.column] <- feature.names[
        na.features,
        replacement.column
      ]
    }
    if (unique.features) {
      fcols <- ncol(x = feature.names)
      if (fcols < gene.column) {
        stop(paste0(
          "gene.column was set to ", gene.column,
          " but feature.tsv.gz (or genes.tsv) only has ",
          fcols, " columns.", " Try setting the gene.column argument to a value <= to ",
          fcols, "."
        ))
      }
      rownames(x = data) <- make.unique(names = feature.names[
        ,
        gene.column
      ])
    }
    if (ncol(x = feature.names) > 2) {
      data_types <- factor(x = feature.names$V3)
      lvls <- levels(x = data_types)
      if (length(x = lvls) > 1 && length(x = full.data) ==
        0) {
        message("10X data contains more than one type and is being returned as a list containing matrices of each type.")
      }
      expr_name <- "Gene Expression"
      if (expr_name %in% lvls) {
        lvls <- c(expr_name, lvls[-which(x = lvls ==
          expr_name)])
      }
      data <- lapply(X = lvls, FUN = function(l) {
        return(data[data_types == l, , drop = FALSE])
      })
      names(x = data) <- lvls
    } else {
      data <- list(data)
    }
    full.data[[length(x = full.data) + 1]] <- data
  }
  list_of_data <- list()
  for (j in 1:length(x = full.data[[1]])) {
    list_of_data[[j]] <- do.call(cbind, lapply(
      X = full.data,
      FUN = `[[`, j
    ))
    list_of_data[[j]] <- as.sparse(x = list_of_data[[j]])
  }
  names(x = list_of_data) <- names(x = full.data[[1]])
  if (length(x = list_of_data) == 1) {
    return(list_of_data[[1]])
  } else {
    return(list_of_data)
  }
}

SeuratObject:::CreateSeuratObject.default
ftype(CreateSeuratObject)
s3_get_method(CreateSeuratObject)
ftype(CreateAssay5Object)
SeuratObject:::CreateAssay5Object


z <- c(1, 3, 4)

pick <- function(i) {
  force(i)
  function(x) x[[i]]
}

pick2 <- function(x, i) {
  x[[i]]
}

pick(3)(z)
x[[1]]

pick2(z, 3)
