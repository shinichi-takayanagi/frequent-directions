# Dummy data
size_col <- 50
size_row <- 10^3
x <- matrix(
  c(rnorm(size_row * size_col), rnorm(size_row * size_col, mean=1)),
  ncol = size_col, byrow = TRUE
)
x <- scale(x)
y <- rep(1:2, each=size_row)
# Show 2D plot using SVD
frequentdirections::plot_svd(x, y)
# Matrix Skethinc(l=6)
b <- frequentdirections::sketching(x, 6, 10^(-8))
# Show 2D plot using sketched matrix and show similar result with the above
# That means that 6 dim is enough to express the original data matrix (x)
frequentdirections::plot_svd(x, y, b)
