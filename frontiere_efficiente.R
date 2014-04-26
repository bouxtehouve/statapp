
frontiere_efficiente<-function (object, frontier = c("both", "lower", "upper"), col = c("black", 
                                                                 "grey"), add = FALSE, labels = TRUE, return = c("mean", "mu"), 
         risk = c("Cov", "Sigma", "CVaR", "VaR"), auto = TRUE, title = TRUE, 
         ...) 
{
  stopifnot(length(col) == 2)
  frontier = match.arg(frontier)
  fullFrontier = frontierPoints(object, frontier = "both", 
                                return = return, risk = risk, auto = auto)
  upperFrontier = frontierPoints(object, frontier = "upper", 
                                 return = return, risk = risk, auto = auto)
  lowerFrontier = frontierPoints(object, frontier = "lower", 
                                 return = return, risk = risk, auto = auto)
  Arg <- match.call(expand.dots = TRUE)
  m <- match(c("xlim", "ylim"), names(Arg), Arg)
  xArg <- as.character(Arg[c(1, m)])[2]
  yArg <- as.character(Arg[c(1, m)])[3]
  if (xArg == "NULL" & yArg == "NULL") {
    yLim = range(fullFrontier[, 2])
    xRange = range(fullFrontier[, 1])
    xDiff = diff(xRange)
    xLim = c(xRange[1] - 2.5 * xDiff/10, xRange[2] + xDiff/10)
    if (!add) {
      if (frontier == "upper" | frontier == "both") {
        plot(upperFrontier, col = col[1], xlim = xLim, 
             ylim = yLim, ann = FALSE, ...)
      }
      else {
        if (frontier == "both") {
          points(fullFrontier, col = col[2], xlim = xLim, 
                 ylim = yLim, ...)
        }
        if (frontier == "lower") {
          plot(lowerFrontier, col = col[2], xlim = xLim, 
               ylim = yLim, ann = FALSE, ...)
        }
      }
    }
    if (frontier == "upper" | frontier == "both") {
      points(upperFrontier, col = col[1], ...)
    }
    if (frontier == "lower" | frontier == "both") {
      points(lowerFrontier, col = col[2], ...)
    }
  }
  else if (xArg != "NULL" & yArg == "NULL") {
    yLim = range(fullFrontier[, 2])
    if (!add) {
      if (frontier == "upper" | frontier == "both") {
        plot(upperFrontier, col = col[1], ylim = yLim, 
             ann = FALSE, ...)
      }
      else {
        if (frontier == "both") {
          points(fullFrontier, col = col[2], ylim = yLim, 
                 ...)
        }
        if (frontier == "lower") {
          plot(fullFrontier, col = col[2], ylim = yLim, 
               ann = FALSE, ...)
        }
      }
    }
    if (frontier == "upper" | frontier == "both") {
      points(upperFrontier, col = col[1], ...)
    }
    if (frontier == "lower" | frontier == "both") {
      points(lowerFrontier, col = col[2], ...)
    }
  }
  else if (xArg == "NULL" & yArg != "NULL") {
    xRange = range(fullFrontier[, 1])
    xDiff = diff(xRange)
    xLim = c(xRange[1] - 2.5 * xDiff/10, xRange[2] + xDiff/10)
    if (!add) {
      if (frontier == "upper" | frontier == "both") {
        plot(upperFrontier, col = col[1], xlim = xLim, 
             ann = FALSE, ...)
      }
      else {
        if (frontier == "both") {
          points(fullFrontier, col = col[2], xlim = xLim, 
                 ...)
        }
        if (frontier == "lower") {
          plot(lowerFrontier, col = col[2], xlim = xLim, 
               ann = FALSE, ...)
        }
      }
    }
    if (frontier == "upper" | frontier == "both") {
      points(upperFrontier, col = col[1], ...)
    }
    if (frontier == "lower" | frontier == "both") {
      points(lowerFrontier, col = col[2], ...)
    }
  }
  else if (xArg != "NULL" & yArg != "NULL") {
    if (!add) {
      if (frontier == "upper" | frontier == "both") {
        plot(fullFrontier, type = "n", ann = FALSE, ...)
        points(upperFrontier, col = col[1], ...)
      }
      if (frontier == "both") {
        points(lowerFrontier, col = col[2], ...)
      }
      if (frontier == "lower") {
        plot(lowerFrontier, col = col[2], ann = FALSE, 
             ...)
      }
    }
    else {
      if (frontier == "upper" | frontier == "both") {
        points(upperFrontier, col = col[1], ...)
      }
      if (frontier == "lower" | frontier == "both") {
        points(lowerFrontier, col = col[2], ...)
      }
    }
  }
  if (labels) {
    mtext(paste(getType(object), "|", getSolver(object)), 
          side = 4, adj = 0, col = "grey", cex = 0.7)
  }
  if (title) {
    labs = attr(fullFrontier, "control")
    title(main = "Frontière efficiente", xlab = paste("Risque[", 
                                                    labs[1], "]", sep = ""), ylab = paste("Rendement[", 
                                                                                          labs[2], "]", sep = ""))
  }
  invisible(fullFrontier)
}

frontiere_efficiente1<-function (object, return = c("mean", "mu"), risk = c("Cov", "Sigma", 
                                                             "CVaR", "VaR"), mText = NULL, col = NULL, xlim = NULL, ylim = NULL, 
                  twoAssets = FALSE) 
{
  offset = 0.1
  risk <- match.arg(risk)
  if (is.null(xlim)) {
    if (risk == "Cov") {
      xmax = max(sqrt(diag(getCov(object))))
    }
    if (risk == "Sigma") {
      xmax = max(sqrt(diag(getSigma(object))))
    }
    if (risk == "CVaR") {
      alpha = getAlpha(object)
      quantiles = colQuantiles(getSeries(object), prob = alpha)
      n.max = which.max(-quantiles)
      r = getSeries(object)[, n.max]
      r = r[r < quantiles[n.max]]
      xmax = -mean(r)
    }
    if (risk == "VaR") {
      xmax = max(-colQuantiles(getSeries(object), prob = alpha))
    }
    xlim = c(0, xmax)
    Xlim = c(xlim[1] - diff(xlim) * offset, xlim[2] + diff(xlim) * 
               offset)
  }
  if (is.null(ylim)) {
    ylim = range(getMean(object))
    Ylim = c(ylim[1] - diff(ylim) * offset, ylim[2] + diff(ylim) * 
               offset)
  }
  frontiere_efficiente(object, return = return, risk = risk, auto = FALSE, 
               xlim = Xlim, ylim = Ylim, pch = 19)
  if (is.null(mText)) 
    mText = getTitle(object)
  mtext(mText, side = 3, line = 0.5, font = 2)
  grid()
  abline(h = 0, col = "grey")
  abline(v = 0, col = "grey")
  data = getData(object)
  spec = getSpec(object)
  constraints = getConstraints(object)
  mvPortfolio = minvariancePortfolio(data, spec, constraints)
  minvariancePoints(object, return = return, risk = risk, auto = FALSE, 
                    pch = 19, col = "red")
  tangencyPoints(object, return = return, risk = risk, auto = FALSE, 
                 pch = 19, col = "blue")
  tangencyLines(object, return = return, risk = risk, auto = FALSE, 
                col = "blue")
  xy = equalWeightsPoints(object, return = return, risk = risk, 
                          auto = FALSE, pch = 15, col = "grey")
  text(xy[, 1] + diff(xlim)/20, xy[, 2] + diff(ylim)/20, "EWP", 
       font = 2, cex = 0.7)
  if (is.null(col)) 
    col = rainbow(6)
  xy = singleAssetPoints(object, return = return, risk = risk, 
                         auto = FALSE, cex = 1.5, col = col, lwd = 2)
  text(xy[, 1] + diff(xlim)/20, xy[, 2] + diff(ylim)/20, rownames(xy), 
       font = 2, cex = 0.7)
  if (twoAssets) {
    twoAssetsLines(object, return = return, risk = risk, 
                   auto = FALSE, lty = 3, col = "grey")
  }
  sharpeRatioLines(object, return = return, risk = risk, auto = FALSE, 
                   col = "orange", lwd = 2)
  invisible(object)
}
