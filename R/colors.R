# Code in this file is adapted from
# https://github.com/btupper/catecolors
# by Ben Tupper.

# K. L. Kelly. Twenty two colors of maximum contrast. Color
# Engineering, 3:26-27, 1965.
# http://www.iscc.org/pdf/PC54_1724_001.pdf
kelly <- function (index, ...) {
  if (missing(index))
    index <- seq_len(nrow(KELLYLUT))
  return(get_lut(KELLYLUT,index,...))
}

# Retrieve one or more of Glasbey et al 256 color specifications.
glasbey <- function (index, ...) {
  if (missing(index))
    index <- seq_len(nrow(GLASBEYLUT))
  return(get_lut(GLASBEYLUT,index,...))
}

# Retrieve one or more of color specifications as hex, rgb triplets or
# a data.frame.
get_lut <- function (LUT, index, 
                     form = c("hex", "rgb", "data.frame")[1],
                     name = FALSE) {
  if (missing(LUT))
    stop("LUT is required")
  if (name)
    nm <- rownames(LUT)
  if (missing(index))
    index <- seq_len(nrows(LUT))
  form <- tolower(form[1])
  if (form == "hex"){
    x <- LUT[index,"hex"]
    if (name)
      names(x) <- rownames(LUT)[index]
  } else if (form == "rgb") {
    x <- as.matrix(LUT[index,c("red","green","blue")])
    if (name)
      names(x) <- rownames(LUT)[index]
  } else
    x <- LUT[index,]
  return(x)
}
