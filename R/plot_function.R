plotting <- function(plot.object, descr.object, factor, col, pch, legendpos, ...){

  nf <- plot.object$nf
  color <- col
  # one-way 
  if(nf == 1){
    xmax <- length(plot.object$levels[[1]])
    y <- descr.object[, 3]
    li <- descr.object[, 5]
    ui <- descr.object[, 6]
    
    plotrix::plotCI(x = 1:xmax, y, li = li,
                    ui = ui, xlim = c(0.8, length(plot.object$levels[[1]]) + 0.3), xaxt = "n", col = color[1], pch = pch[1], ...)
    axis(side = 1, at = 1:1:length(plot.object$levels[[1]]), labels = plot.object$levels[[1]], ...)
  } else {
    
    TYPE <- plot.object$Type    
    Faktor <- strsplit(factor, ":")[[1]]
    nadat2 <- plot.object$nadat2
    levels <- plot.object$levels
    mu <- plot.object$mu
    lower <- plot.object$lower
    upper <- plot.object$upper
    lev_names <- plot.object$lev_names
    mean_out <- descr.object$Means
    CI <- matrix(c(descr.object$Lower, descr.object$Upper), ncol = 2)
    colnames(CI) <- c("CIl", "CIu")
    fac_names <- plot.object$fac_names
    fl <- plot.object$fl
    fac_names_original <- plot.object$fac_names_original

    # plots for interactions
    if (TYPE == "nested") {
      # main effect
      if (factor == nadat2[1]) {
        plotrix::plotCI(x = 1:length(levels[[1]]), mu[[1]],
                        li = lower[[1]], ui = upper[[1]], xlim = c(0.8, length(levels[[1]]) + 0.3),
                        col = color[1], pch = pch[1], xaxt = "n", ...)
        axis(side = 1, at = 1:1:length(levels[[1]]), labels = levels[[1]], ...)
      } else if (factor == fac_names[2] && nf == 2) {
        plotrix::plotCI(x = 1:length(levels[[2]]), mean_out,
                        li = CI[, 1], ui = CI[, 2], xlim = c(0.8, length(levels[[2]]) + 0.3),
                        ylim = c(min(CI) - 1, max(CI) + 1), col = color[1], pch = pch[1], xaxt = "n", ...)
        axis(side = 1, at = 1:1:length(levels[[2]]), labels = lev_names[, 2], ...)
        aa <- length(levels[[2]]) / fl[1] ^ 2
        bb <- length(levels[[2]]) / fl[1]
        cc <- length(levels[[2]]) + aa
        ss <- seq(from = - aa, to = cc, by = bb)
        axis(side = 3, at = ss[2:(length(ss) - 1)], labels = levels[[1]], ...)
      } else if(factor %in% fac_names && nf == 3 && !(factor == nadat2[1])) {
        error <- "For three-way nested design, only the main effect
    can be plotted."
      }}
    
    
    if (TYPE == "crossed") {
      # plot of main effects
      for (i in 1:nf) {
        if (factor == nadat2[i]) {
          plotrix::plotCI(x = 1:length(levels[[i]]), mu[[i]],
                          li = lower[[i]], ui = upper[[i]], xlim = c(0.8, length(levels[[i]]) + 0.3),
                          col = color[1], pch = pch[1], xaxt = "n", ...)
          axis(side = 1, at = 1:1:length(levels[[i]]), labels = levels[[i]], ...)
        }
      }
      
      # two-way plots
      fac_names_twofold <- plot.object$fac_names_original[ - (1:nf)]
      fac_names_twofold <- fac_names_twofold[1:choose(nf, 2)]
      
      if (factor %in% fac_names_twofold) {
        dat2 <- plot.object$dat2
        fl <- plot.object$fl
        alpha <- plot.object$alpha
        
        fak1 <- Faktor[1]
        fak2 <- Faktor[2]
        posi <- which(fac_names_original[1:nf] == fak1)
        posi2 <- which(fac_names_original[1:nf] == fak2)
        
        nmu <- matrix(by(dat2[, 1], dat2[, c(fak1, fak2)], mean),
                      nrow = fl[posi])
        nsigma <- matrix(by(dat2[, 1], dat2[, c(fak1, fak2)], var),
                         nrow = fl[posi])
        nn_groups <- matrix(by(dat2[, 1], dat2[, c(fak1, fak2)],
                               length), nrow = fl[posi])
        nlower <- nmu - sqrt(nsigma/ nn_groups) *
          qt(1 - alpha / 2, df = nn_groups)
        nupper <- nmu + sqrt(nsigma / nn_groups) *
          qt(1 - alpha / 2, df = nn_groups)
        
        
        plotrix::plotCI(x = 1:dim(nmu)[2],
                        nmu[1, ],
                        li = nlower[1, ],
                        ui = nupper[1, ], xlim = c(0.8, (dim(nmu)[2] + 0.3)),
                        ylim = c(min(nlower) - 1, max(nupper) + 1),
                        col = color[1], pch = pch[1], xaxt = "n", ...)
        axis(side = 1, at = (1:1:dim(nmu)[2]),
             labels = unlist(levels[posi2]), ...)
        for (i in 2:dim(nmu)[1]) {
          plotrix::plotCI(x = ((1:dim(nmu)[2]) + 0.07 * i),
                          nmu[i, ], li = nlower[i, ],
                          ui = nupper[i, ], add = TRUE, col = color[i], pch = pch[1], ...)
          
        }
        legend(legendpos, 
               legend = unlist(levels[posi]),
               col = color[1:dim(nmu)[1]],
               seg.len = 0.5, 
               lty = rep(1, dim(nmu)[1]))
        
      } else if (length(Faktor) == 3) {
        # three-way plots
        
        fak1 <- Faktor[1]
        fak2 <- Faktor[2]
        fak3 <- Faktor[3]
        
        posi1 <- which(fac_names_original[1:nf] == fak1)
        posi2 <- which(fac_names_original[1:nf] == fak2)
        posi3 <- which(fac_names_original[1:nf] == fak3)
        
        if (nf == 3){          
          mu3 <- descr.object$Means
          lower3 <- descr.object[, 7]
          upper3 <- descr.object[, 8]
        } else {
          # c(t(...)) to get correct order of factors, analogous to case above
          mu3 <- c(t(matrix(by(dat2[, 1], dat2[, c(fak1, fak3, fak2)], mean),
                            nrow = fl[posi1])))
          nsigma <- c(t(matrix(by(dat2[, 1], dat2[, c(fak1, fak3, fak2)], var),
                               nrow = fl[posi1])))
          nn_groups <- c(t(matrix(by(dat2[, 1], dat2[, c(fak1, fak3, fak2)],
                                     length), nrow = fl[posi1])))
          lower3 <- mu3 - sqrt(nsigma/ nn_groups) *
            qt(1 - alpha / 2, df = nn_groups)
          upper3 <- mu3 + sqrt(nsigma / nn_groups) *
            qt(1 - alpha / 2, df = nn_groups)          
        }
        
        delta <- seq(from = 0, by = 0.05,
                     length = (fl[posi1] * fl[posi2] + 1))
        
        plotrix::plotCI(x = 1:fl[posi3],
                        mu3[1:fl[posi3]],
                        li = lower3[1:fl[posi3]],
                        ui = upper3[1:fl[posi3]], xlim = c(0.8, (fl[posi3] + 0.3)),
                        ylim = c(min(lower3) - 1, max(upper3) + 1),
                        col = color[1], pch = pch[1], xaxt = "n", ...)
        axis(side = 1, at = 1:1:fl[posi3], labels = levels[[posi3]], ...)
        
        color2 <- rep(color[1:fl[posi2]], fl[posi1])
        pch2 <- rep(pch[1:fl[posi1]], each = fl[posi2])
        
        for (j in 1:(fl[posi1]*fl[posi2]-1)) {
          start <- fl[posi3]*j+1
          stopp <- fl[posi3]*(j+1)
          plotrix::plotCI(x = ((1:fl[posi3]) + delta[j+1]),
                          mu3[start:stopp], li = lower3[start:stopp],
                          ui = upper3[start:stopp], add = TRUE, col = color2[j+1], pch = pch2[j+1], ...)
        }
        legend(legendpos, 
               legend = c(fac_names_original[[posi1]], levels[[posi1]], fac_names_original[[posi2]], levels[[posi2]]), box.lty = 0,
               col = c(0, rep(1, fl[posi1]), 0, color[1:fl[posi2]]),
               pch = c(NA, pch[1:fl[posi1]], NA, rep(NA, fl[posi2])),
               seg.len = 0.5,
               lty = c(NA, rep(NA, fl[posi1]), NA, rep(1, fl[posi2])))
        
      } else if (length(Faktor) >= 4) {
        stop("Higher-way interactions cannot be plotted!")
      }
    }
  }
}
