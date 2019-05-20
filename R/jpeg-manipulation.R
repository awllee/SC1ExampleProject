library(jpeg)

# input: a 3D array, essentially 3 matrices of equal size
# output: a matrix of grey-scale values
grayscale <- function(img) {
  dims <- dim(img)
  stopifnot(length(dims) == 3 && dims[3] == 3)
  mtx1 <- img[,,1]
  mtx2 <- img[,,2]
  mtx3 <- img[,,3]
  mtx <- (mtx1 + mtx2 + mtx3)/3
  return(mtx)
}

# input: name of a JPEG file
# output: a raster array for the image
load.jpeg <- function(filename) {
  img <- readJPEG(filename)
  return(img)
}

# input: raster array for the image
# effect: image is plotted
view.image <- function(img) {
  rst <- as.raster(img)
  plot(rst)
}

# input: matrix
# output: matrix with entries thresholded below/above at 0/1
fix.image <- function(mtx) {
  return(pmin(pmax(mtx,0),1))
}
