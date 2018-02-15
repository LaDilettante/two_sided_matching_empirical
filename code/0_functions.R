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

f_plot_mcmc <- function(x, main = "") {
  plot(x, type = "l", main = main)
  lines(1:length(x), zoo::rollmean(x, k = 50, fill = NA), col = "red")
  plot(density(x), main = main)
}

my.plot.mcmc <- function (x, trace = TRUE, density = TRUE, smooth = FALSE, bwf, 
                          auto.layout = TRUE, ask = FALSE, parameters, ...) 
{
  oldpar <- NULL
  on.exit(par(oldpar))
  if (auto.layout) {
    mfrow <- coda:::set.mfrow(Nchains = nchain(x), Nparms = nvar(x), 
                              nplots = trace + density)
    oldpar <- par(mfrow = mfrow)
  }
  for (i in 1:nvar(x)) {
    y <- mcmc(as.matrix(x)[, i, drop = FALSE], start(x), 
              end(x), thin(x))
    if (trace) 
      traceplot(y, smooth = smooth, ...) ; abline(h=parameters[i], col="red")
    if (density) {
      if (missing(bwf)) {
        densplot(y, ...); abline(v=parameters[i], col="red")
      } else densplot(y, bwf = bwf, ...)
    }
    if (i == 1) 
      oldpar <- c(oldpar, par(ask = ask))
  }
}
