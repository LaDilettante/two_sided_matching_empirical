# Install and load packages
f_install_and_load <- function(packs) {
  new.packs <- packs[!(packs %in% installed.packages()[ ,"Package"])]
  lapply(new.packs, install.packages, repos="http://cran.rstudio.com/", dependencies=TRUE)
  lapply(packs, library, character.only=TRUE)
}

heatmap2 <- function(x) {
  heatmap(x, scale = "none", Rowv = NA, Colv = NA, revC = T)
}

postmean <- function(x, start = nrow(x) / 2) {
  colMeans(x[start:nrow(x), ])
}

postsd <- function(x, start = nrow(x) / 2) {
  apply(x[start:nrow(x), ], 1, sd)
}