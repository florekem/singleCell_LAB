library(sloop)
library(lobstr)

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
