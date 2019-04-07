all_zero_row_index <- function(x, eps){
  which(abs(rowSums(x)) < eps)
}

#' Compute a sketch matrix of input matrix
#'
#' Compute a sketch matrix of input matrix
#'
#' @param a Original matrix to be sketched (n x m)
#' @param l The number of rows in sketched matrix (l x m)
#' @param eps If a value is smaller than eps, that is considered as equal to zero. The default value is 10^(-8)
#' @export
sketching <- function(a, l, eps=10^(-8)){
  m <- ncol(a)
  n <- nrow(a)
  # Input error handling
  if(floor(l / 2) >= m){stop("l must be smaller than m * 2")}
  if(l >= n){stop("l must not be greater than n")}

  b <- matrix(0, nrow = l, ncol = m)
  zero_row_index <- all_zero_row_index(b, eps)
  for(i in seq_len(l)){
    # Fill first all zero row by a[i,]
    b[zero_row_index[1], ] <- a[i, ]
    #Remove first element because we already used it
    zero_row_index <- tail(zero_row_index, -1)
    if(length(zero_row_index) == 0){
      b_svd <- svd(b)
      v <- b_svd$v
      sigma <- b_svd$d
      delta <- sigma[floor(l/2)]^2
      sigma_tilde <- sqrt(pmax(sigma^2 - delta, 0))
      b <- diag(sigma_tilde) %*% t(v)
      # Update zero lists
      zero_row_index <- all_zero_row_index(b, eps)
    }
  }
  b
}
#' Plot data
#'
#' Plot data
#'
#' @param data x
#' @param label y
#' @param x skeched matrix
#' @export
plot_svd <- function(data, label, x = data){
  v <- svd(x)$v
  # Not cool code...
  if(sum(v[,1]) <= 0){v[,1] <- -v[,1]}
  if(sum(v[,2]) <= 0){v[,2] <- -v[,2]}
  # Projection matrix(x_p = XV = UΣ) and plot
  x_p <- data %*% v[,1:2]
  plot(x_p[,1], x_p[,2], col=label, pch=16)
}
