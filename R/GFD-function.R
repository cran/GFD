par_old <- par()$oma
on.exit(par(oma = par_old))

#' Tests for general factorial designs
#' 
#' The GFD function calculates the Wald-type statistic (WTS), the ANOVA-type 
#' statistic (ATS) as well as a permutation version of the WTS for general 
#' factorial designs.
#' 
#' @param formula A model \code{\link{formula}} object. The left hand side
#'   contains the response variable and the right hand side contains the factor
#'   variables of interest. An interaction term must be specified.
#' @param data A data.frame, list or environment containing the variables in 
#'   \code{formula}. The default option is \code{NULL}.
#' @param nperm The number of permutations used for calculating the permuted 
#'   Wald-type statistic. The default option is 10000.
#' @param alpha A number specifying the significance level; the default is 0.05.
#' @param plot_CI An indicator for whether a plot of the results should be 
#'   shown; the default is \code{FALSE}.
#' @param Title A title for the plot. The default is \code{NULL}.
#' @param line_width A number specifying the \code{lwd}-option for the plots.
#'   The default is 2.
#'   
#' @details The package provides the Wald-type statistic, a permuted version
#'   thereof as well as the ANOVA-type statistic for general factorial designs,
#'   even with non-normal error terms and/or heteroscedastic variances. It is
#'   implemented for both crossed and hierarchically nested designs and allows
#'   for an arbitrary number of factor combinations as well as different sample
#'   sizes in the crossed design.

#'   
#' @return A \code{GFD} object containing the following components:
#' \item{Descriptive}{Some descriptive statistics of the data for all factor
#'   level combinations. Displayed are the number of individuals per factor
#'   level combination, the mean, variance and 100*(1-alpha)\% confidence
#'   intervals.}
#'  \item{WTS}{The value of the WTS along with degrees of freedom of the central chi-square distribution and p-value, as well as the p-value of the permutation procedure.}
#'  \item{ATS}{The value of the ATS, degrees of freedom of the central F distribution and the corresponding p-value.}
#' 
#' @examples
#' GFD(weightgain ~ source * type, data = HSAUR::weightgain)
#' 
#' data(startup)
#' model <- GFD(Costs ~ company, data = startup)
#' summary(model)
#' 
#' @references Friedrich, S., Konietschke, F., Pauly, M.(2015). GFD - An R-package
#' for the Analysis of General Factorial Designs - along with a Graphical User Interface. Submitted to Journal of Statistical Software.
#' 
#' Pauly, M., Brunner, E., Konietschke, F.(2015). Asymptotic Permutation Tests in General Factorial Designs. Journal of the Royal Statistical Society - Series B 77, 461-473.
#' 
#' @import RGtk2
#' 
#' @importFrom graphics axis legend par plot title
#' @importFrom stats ecdf formula model.frame pchisq pf qt terms var
#' @importFrom utils read.table
#' 
#' @export

GFD <- function(formula, data = NULL, nperm = 10000,
                alpha = 0.05, plot_CI = FALSE, Title = NULL,
                line_width = 2){
  
  input_list <- list(formula = formula, data = data,
                     nperm = nperm, alpha = alpha,
                     plot_CI = plot_CI, Title = Title, line_width = line_width)
  dat <- model.frame(formula, data)
  subject <- 1:nrow(dat)
  dat2 <- data.frame(dat, subject = subject)
  nf <- ncol(dat) - 1
  nadat <- names(dat)
  nadat2 <- nadat[-1]
  fl <- NA
  for (aa in 1:nf) {
    fl[aa] <- nlevels(as.factor(dat[ ,aa + 1]))
  }
  levels <- list()
  for (jj in 1:nf) {
    levels[[jj]] <- levels(as.factor(dat[ ,jj + 1]))
  }
  lev_names <- expand.grid(levels)
  
  if (nf == 1) {
    # one-way layout
    dat2 <- dat2[order(dat2[, 2]), ]
    response <- dat2[, 1]
    nr_hypo <- attr(terms(formula), "factors")
    fac_names <- colnames(nr_hypo)
    n <- plyr::ddply(dat2, nadat2, plyr::summarise, Measure = length(subject),
                     .drop = F)$Measure
    # contrast matrix
    hypo <- diag(fl) - matrix(1 / fl, ncol = fl, nrow = fl)
    WTS_out <- matrix(NA, ncol = 3, nrow = 1)
    ATS_out <- matrix(NA, ncol = 4, nrow = 1)
    WTPS_out <- rep(NA, 1)
    rownames(WTS_out) <- fac_names
    rownames(ATS_out) <- fac_names
    names(WTPS_out) <- fac_names
    results <- Stat(data = response, n = n, hypo, nperm = nperm, alpha)
    WTS_out <- results$WTS
    ATS_out <- results$ATS
    WTPS_out <- results$WTPS
    mean_out <- results$Mean
    Var_out <- results$Cov
    CI <- results$CI
    colnames(CI) <- c("CIl", "CIu")
    descriptive <- cbind(lev_names, n, mean_out, Var_out, CI)
    colnames(descriptive) <- c(nadat2, "n", "Means", "Variances",
                               paste("Lower", 100 * (1 - alpha), "%", "CI"),
                               paste("Upper", 100 * (1 - alpha), "%", "CI"))
    
    if (plot_CI == TRUE) {
      GUIplotOneWay <- function() {
        plotting <- function(button, user.data) {
          error <- NULL
          Title <- filename2$getText()
          line_width <- as.numeric(filename3$getText())
          if (!is.null(error)) {
            hbox <- RGtk2::gtkHBoxNew()
            vbox$packStart(hbox, FALSE, FALSE, 0)
            label <- RGtk2::gtkLabel(error)
            hbox$packStart(label, FALSE, FALSE, 0)
          }
          plotrix::plotCI(x = 1:length(levels[[1]]), mean_out, li = descriptive[, 5],
                          ui = descriptive[, 6], xlab = fac_names, ylab = "Means",
                          col = 1, lwd = line_width, xaxt = "n",
                          xlim = c(0.8, length(levels[[1]]) + 0.3), main = Title,
                          cex.axis = 1.3, cex.lab = 1.3, font.axis = 2, font.lab = 2)
          axis(side = 1, at = 1:1:length(levels[[1]]), labels = levels[[1]], las = 0, cex.axis = 1.3, cex.lab = 1.3, font.axis = 2)
        }
        # Create window
        window <- RGtk2::gtkWindow()
        # Add title
        window["title"] <- "Plot"
        # Add a frame
        frame <- RGtk2::gtkFrameNew("Please choose the parameters for your plot.")
        window$add(frame)
        # Create vertical container for file name entry
        vbox <- RGtk2::gtkVBoxNew(FALSE, 8)
        vbox$setBorderWidth(24)
        frame$add(vbox)
        # Add horizontal container for every widget line
        hbox <- RGtk2::gtkHBoxNew(FALSE, 8)
        vbox$packStart(hbox, FALSE, FALSE, 0)
        # Add an horizontal container to specify parameters
        hbox <- RGtk2::gtkHBoxNew(FALSE, 8)
        vbox$packStart(hbox, FALSE, FALSE, 0)
        label2 <- RGtk2::gtkLabelNewWithMnemonic("_Title")
        hbox$packStart(label2, FALSE, FALSE, 0)
        # Add entry in the second column; named "filename2"
        filename2 <- RGtk2::gtkEntryNew()
        filename2$setWidthChars(10)
        label2$setMnemonicWidget(filename2)
        hbox$packStart(filename2, FALSE, FALSE, 0)
        label3 <- RGtk2::gtkLabelNewWithMnemonic("_lwd")
        hbox$packStart(label3, FALSE, FALSE, 0)
        # Add entry in the second column; named "filename3"
        filename3 <- RGtk2::gtkEntryNew()
        filename3$setWidthChars(10)
        filename3$setText(2)
        label3$setMnemonicWidget(filename3)
        hbox$packStart(filename3, FALSE, FALSE, 0)
        # Add button
        the.buttons <- RGtk2::gtkHButtonBoxNew()
        the.buttons$setBorderWidth(5)
        vbox$add(the.buttons)
        the.buttons$setLayout("spread")
        the.buttons$setSpacing(40)
        buttonOK <- RGtk2::gtkButtonNewFromStock("gtk-ok")
        RGtk2::gSignalConnect(buttonOK, "clicked", plotting)
        the.buttons$packStart(buttonOK,fill=F)
        buttonCancel <- RGtk2::gtkButtonNewFromStock("gtk-close")
        RGtk2::gSignalConnect(buttonCancel, "clicked", window$destroy)
        the.buttons$packStart(buttonCancel, fill=F)
      }
      GUIplotOneWay()
    }
    WTS_output <- c(WTS_out, WTPS_out)
    names(WTS_output) <- cbind ("Test statistic", "df",
                                "p-value", "p-value WTPS")
    names(ATS_out) <- cbind("Test statistic", "df1", "df2", "p-value")
    output <- list()
    output$input <- input_list
    output$Descriptive <- descriptive
    output$WTS <- WTS_output
    output$ATS <- ATS_out
    # end one-way layout ------------------------------------------------------
  } else {
    lev_names <- lev_names[do.call(order, lev_names[, 1:nf]), ]
    # sorting data according to factors
    dat2 <- dat2[do.call(order, dat2[, 2:(nf + 1)]), ]
    response <- dat2[, 1]
    nr_hypo <- attr(terms(formula), "factors")
    fac_names <- colnames(nr_hypo)
    fac_names_original <- fac_names
    perm_names <- t(attr(terms(formula), "factors")[-1, ])
    n <- plyr::ddply(dat2, nadat2, plyr::summarise, Measure = length(subject),
                     .drop = F)$Measure
    if (length(fac_names) == nf) {
      # delete factorcombinations which don't exist
      n <- n[n != 0]
    }
    # mixture of nested and crossed designs is not possible
    if (length(fac_names) != nf && 2 %in% nr_hypo) {
      stop("A model involving both nested and crossed factors is
           not impemented!")
    }
    # only 3-way nested designs are possible
    if (length(fac_names) == nf && nf >= 4) {
      stop("Four- and higher way nested designs are
           not implemented!")
    }
    # no factor combinations with less than 2 observations
    if (0 %in% n || 1 %in% n) {
      stop("There is at least one factor-level combination
           with less than 2 observations!")
    }
    # correct labeling of factors in nested design
    if (length(fac_names) == nf) {
      if (nf == 2) {
        if (all(levels(as.factor(dat2[, 3][dat2[, 2] == levels[[1]][1]]))
                == levels(as.factor(dat2[, 3][dat2[, 2] == levels[[1]][2]])))) {
          stop("The levels of the nested factor must be
               named without repetitions!")
        }
      } else if (nf == 3) {
        if (all(levels(as.factor(dat2[, 3][dat2[, 2] == levels[[1]][1]]))
                == levels(as.factor(dat2[, 3][dat2[, 2] == levels[[1]][2]]))) ||
              all(levels(as.factor(dat2[, 4][dat2[, 3] == levels[[2]][1]]))
                  == levels(as.factor(dat2[, 4][dat2[, 3] == levels[[2]][fl[2] / fl[1] + 1]])))) {
          stop("The levels of the nested factor must be
                 named without repetitions!")
        }
      }
    }
    if (length(fac_names) == nf) {
      # nested
      TYPE <- "nested"
      hypo_matrices <- HN(fl)
      # create correct level combinations
      blev <- list()
      lev_names <- list()
      for (ii in 1:length(levels[[1]])) {
        blev[[ii]] <- levels(as.factor(dat[, 3][dat[, 2] == levels[[1]][ii]]))
        lev_names[[ii]] <- rep(levels[[1]][ii], length(blev[[ii]]))
      }
      if (nf == 2) {
        lev_names <- as.factor(unlist(lev_names))
        blev <- as.factor(unlist(blev))
        lev_names <- cbind.data.frame(lev_names, blev)
      } else {
        lev_names <- lapply(lev_names, rep,
                            length(levels[[3]]) / length(levels[[2]]))
        lev_names <- lapply(lev_names, sort)
        lev_names <- as.factor(unlist(lev_names))
        blev <- lapply(blev, rep, length(levels[[3]]) / length(levels[[2]]))
        blev <- lapply(blev, sort)
        blev <- as.factor(unlist(blev))
        lev_names <- cbind.data.frame(lev_names, blev, as.factor(levels[[3]]))
      }
    } else {
      # crossed
      TYPE <- "crossed"
      hypo_matrices <- HC(fl, perm_names, fac_names)[[1]]
      fac_names <- HC(fl, perm_names, fac_names)[[2]]
    }
    if (length(fac_names) != length(hypo_matrices)) {
      stop("Something is wrong: Perhaps a missing interaction term in formula?")
    }
    WTS_out <- matrix(NA, ncol = 3, nrow = length(hypo_matrices))
    ATS_out <- matrix(NA, ncol = 4, nrow = length(hypo_matrices))
    WTPS_out <- rep(NA, length(hypo_matrices))
    rownames(WTS_out) <- fac_names
    rownames(ATS_out) <- fac_names
    names(WTPS_out) <- fac_names
    colnames(ATS_out) <- c("Test statistic", "df1", "df2", "p-value")
    # calculate results
    for (i in 1:length(hypo_matrices)) {
      results <- Stat(data = response, n = n, hypo_matrices[[i]],
                      nperm = nperm, alpha)
      WTS_out[i, ] <- results$WTS
      ATS_out[i, ] <- results$ATS
      WTPS_out[i] <- results$WTPS
    }
    mean_out <- results$Mean
    Var_out <- results$Cov
    CI <- results$CI
    colnames(CI) <- c("CIl", "CIu")
    descriptive <- cbind(lev_names, n, mean_out, Var_out, CI)
    colnames(descriptive) <- c(nadat2, "n", "Means", "Variances",
                               paste("Lower", 100 * (1 - alpha),"%", "CI"),
                               paste("Upper", 100 * (1 - alpha),"%", "CI"))
    
    # calculate group means, variances and CIs ----------------------------
    mu <- list()
    sigma <- list()
    n_groups <- list()
    lower <- list()
    upper <- list()
    for (i in 1:nf) {
      mu[[i]] <- c(by(dat2[, 1], dat2[, i + 1], mean))
      sigma[[i]] <- c(by(dat2[, 1], dat2[, i + 1], var))
      n_groups[[i]] <- c(by(dat2[, 1], dat2[, i + 1], length))
      lower[[i]] <- mu[[i]] - sqrt(sigma[[i]] / n_groups[[i]]) *
        qt(1 - alpha / 2, df = n_groups[[i]])
      upper[[i]] <- mu[[i]] + sqrt(sigma[[i]] / n_groups[[i]]) *
        qt(1 - alpha / 2, df = n_groups[[i]])
    }
    
    # Plotting ----------------------------------------------------
    if (plot_CI == TRUE) {
      calculateGUIplot <- function() {
        plotting <- function(button, user.data) {
          error <- NULL
          error1 <- NULL
          Faktor <- filename$getText()
          Title <- filename2$getText()
          line_width <- as.numeric(filename3$getText())
          
          if (!(Faktor %in% fac_names)) {
            error1 <- "Please enter a valid factor name"
          }
                  
          # plots for interactions
          if (TYPE == "nested") {
            # main effect
            if (Faktor == nadat2[1]) {
              plotrix::plotCI(x = 1:length(levels[[1]]), mu[[1]],
                              li = lower[[1]], ui = upper[[1]],
                              xlab = nadat2[[1]], ylab = "Means",
                              col = 1, lwd = line_width, xaxt = "n",
                              xlim = c(0.8, length(levels[[1]]) + 0.3),
                              main = Title, font.lab = 2,
                              cex.axis = 1.3, cex.lab = 1.3, font.axis = 2)
              axis(side = 1, at = 1:1:length(levels[[1]]), labels = levels[[1]],
                   las = 0, cex.axis = 1.3, cex.lab = 1.3, font.axis = 2)
            } else if (Faktor == fac_names[2] && nf == 2) {
              for_plots <- cbind(lev_names, mean_out, CI)
              plotrix::plotCI(x = 1:length(levels[[2]]), for_plots[, 3],
                              li = for_plots[, 4], ui = for_plots[, 5],
                              xlab = fac_names[2], ylab = "Means",
                              col = 1, lwd = line_width,
                              ylim = c(min(CI) - 1, max(CI) + 1),
                              xaxt = "n", xlim = c(0.8, length(levels[[2]]) + 0.3),
                              cex.axis = 1.3, cex.lab = 1.3,
                              font.axis = 2, font.lab = 2)
              title(Title, line = 3)
              axis(side = 1, at = 1:1:length(levels[[2]]), labels = lev_names[, 2],
                   las = 0, cex.axis = 1.3, cex.lab = 1.3, font.axis = 2)
              aa <- length(levels[[2]]) / fl[1] ^ 2
              bb <- length(levels[[2]]) / fl[1]
              cc <- length(levels[[2]]) + aa
              ss <- seq(from = - aa, to = cc, by = bb)
              axis(side = 3, at = ss[2:(length(ss) - 1)], labels = levels[[1]],
                   las = 1, cex.axis = 1.3, cex.lab = 1.3, font.axis = 2)
            } else if(Faktor %in% fac_names && nf == 3 && !(Faktor == nadat2[1])) {
              error <- "For three-way nested design, only the main effect
                      can be plotted."
            }}
          if (TYPE == "crossed") {
            # plot of main effects
            for (i in 1:nf) {
              if (Faktor == nadat2[i]) {
                plotrix::plotCI(x = 1:length(levels[[i]]), mu[[i]],
                                li = lower[[i]], ui = upper[[i]],
                                xlab = nadat2[[i]], ylab = "Means",
                                col = 1, lwd = line_width, xaxt = "n",
                                xlim = c(0.8, length(levels[[i]]) + 0.3),
                                main = Title, font.lab = 2,
                                cex.axis = 1.3, cex.lab = 1.3, font.axis = 2)
                axis(side = 1, at = 1:1:length(levels[[i]]), labels = levels[[i]],
                     las = 0, cex.axis = 1.3, cex.lab = 1.3, font.axis = 2)
              }}

            # two-fold interactions for three- and higher-way layout
            fac_names_twofold <- fac_names_original[ - (1:nf)]
            fac_names_twofold <- fac_names_twofold[1:choose(nf, 2)]
            
            if (Faktor %in% fac_names_twofold) {
              nmu <- list()
              nsigma <- list()
              nn_groups <- list()
              nupper <- list()
              nlower <- list()
              new_levels <- list()
              counter <- 1
              for (i in 2:nf) {
                for (j in (i + 1):(nf + 1)) {
                  nmu[[counter]] <- matrix(by(dat2[, 1], dat2[, c(i, j)], mean),
                                           nrow = fl[i - 1])
                  nsigma[[counter]] <- matrix(by(dat2[, 1], dat2[, c(i, j)], var),
                                              nrow = fl[i - 1])
                  nn_groups[[counter]] <- matrix(by(dat2[, 1], dat2[, c(i, j)],
                                                    length), nrow = fl[i - 1])
                  nlower[[counter]] <- nmu[[counter]] -
                    sqrt(nsigma[[counter]] / nn_groups[[counter]]) *
                    qt(1 - alpha / 2, df = nn_groups[[counter]])
                  nupper[[counter]] <- nmu[[counter]] +
                    sqrt(nsigma[[counter]] / nn_groups[[counter]]) *
                    qt(1 - alpha / 2, df = nn_groups[[counter]])
                  new_levels[[counter]] <- list(levels[[i - 1]], levels[[j - 1]])
                  counter <- counter + 1
                }
              }
              names(nmu) <- fac_names_twofold
              names(nupper) <- fac_names_twofold
              names(nlower) <- fac_names_twofold
              place <- which(Faktor == fac_names_twofold)
              xxx <- rep(NA, nf)
              for (ii in 1:nf) {
                xxx[ii] <- grepl(fac_names_original[ii], Faktor, fixed = TRUE)
              }
              pos <- which(xxx == TRUE)[2]
              pos2 <- which(xxx == TRUE)[1]
              par(oma = c(6, 3, 2.5, 2.5))
              plotrix::plotCI(x = 1:length(new_levels[[place]][[2]]),
                              nmu[[Faktor]][1, ],
                              li = nlower[[Faktor]][1, ],
                              ui = nupper[[Faktor]][1, ],
                              xlab = fac_names_original[pos],
                              ylab = "Means", col = 1,
                              lwd = line_width,
                              ylim = c(min(nlower[[Faktor]]) - 1, max(nupper[[Faktor]]) + 1),
                              xaxt = "n",
                              xlim = c(0.8, length(new_levels[[place]][[2]]) + 0.3),
                              main = Title,
                              cex.axis = 1.3, cex.lab = 1.3, font.axis = 2, font.lab = 2)
              axis(side = 1, at = 1:1:length(new_levels[[place]][[2]]),
                   labels = new_levels[[place]][[2]], las = 0, cex.axis = 1.3, cex.lab = 1.3,
                   font.axis = 2)
              for (i in 2:length(new_levels[[place]][[1]])) {
                plotrix::plotCI(x = ((1:length(new_levels[[place]][[2]])) + 0.07 * i),
                                nmu[[Faktor]][i, ], li = nlower[[Faktor]][i, ],
                                ui = nupper[[Faktor]][i, ], add = TRUE, col = i,
                                lwd = line_width, cex.axis = 1.3, cex.lab = 1.3,
                                font.axis = 2, font.lab = 2)
              }
              par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0),
                  new = TRUE)
              plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
              legend("bottom", xpd = TRUE, horiz = TRUE, inset = c(0, 0),
                     box.lty = 0,
                     legend = c(fac_names_original[pos2], new_levels[[place]][[1]]),
                     col = c(1, 1:length(new_levels[[place]][[1]])), lwd = line_width,
                     seg.len = 0.5, text.font = 2,
                     lty = c(NA, rep(1, length(new_levels[[place]][[1]]))))
            } else if (nf == 3 && Faktor == fac_names[length(fac_names)]) {
              # three-way
              Var1 = Var2 = Var3 = CIl = CIu = NULL
              for_plots <- cbind(lev_names, mean_out, CI)
              group <- list()
              for (i in 1:length(levels[[1]])) {
                group[[i]] <- subset(for_plots, Var1 == as.factor(levels[[1]])[i],
                                     select = c(Var2, Var3, mean_out, CIl, CIu))
              }
              next_group <- list()
              new_group <- list()
              for (j in 1:length(group)) {
                for (l in 1:length(levels[[2]])) {
                  next_group[[l]] <- subset(group[[j]],
                                            Var2 == as.factor(levels[[2]])[l],
                                            select = c(mean_out, CIl, CIu))
                }
                new_group[[j]] <- next_group
              }
              counter <- 1
              delta <- seq(from = 0, by = 0.05,
                           length = length(levels[[1]]) * length(levels[[2]]) + 1)
              par(oma = c(6, 3, 2.5, 2.5))
              plotrix::plotCI(x = 1:length(levels[[3]]), new_group[[1]][[1]][, 1],
                              li = new_group[[1]][[1]][, 2],
                              ui = new_group[[1]][[1]][, 3],
                              xlab = fac_names_original[[3]],
                              ylab = "Means", col = 1,
                              lwd = line_width, pch = 2,
                              ylim = c(min(CI) - 1, max(CI) + 1),
                              xaxt = "n", xlim = c(0.8, length(levels[[3]]) + 0.3),
                              main = Title, cex.axis = 1.3, cex.lab = 1.3,
                              font.axis = 2, font.lab = 2)
              axis(side = 1, at = 1:1:length(levels[[3]]), labels = levels[[3]],
                   las = 0, cex.axis = 1.3, cex.lab = 1.3, font.axis = 2)
              for (j in 1:length(levels[[1]])) {
                for (i in 1:length(levels[[2]])) {
                  plotrix::plotCI(x = ((1:length(levels[[3]])) + delta[counter]),
                                  new_group[[j]][[i]][, 1], li = new_group[[j]][[i]][, 2],
                                  ui = new_group[[j]][[i]][, 3], add = TRUE, col = i,
                                  lwd = line_width, pch = 2 * j, cex.axis = 1.3, cex.lab = 1.3,
                                  font.axis = 2, font.lab = 2)
                  counter <- counter + 1
                }}
              par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
              plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
              legend("bottomright", xpd = TRUE, horiz = TRUE, inset = c(0, 0), box.lty = 0,
                     legend = c(levels[[1]], levels[[2]]),
                     col = c(rep(1, length(levels[[1]])), 1:length(levels[[2]])),
                     pch = c(2 * (1:length(levels[[1]])), rep(NA, length(levels[[2]]))),
                     lwd = line_width, seg.len = 0.5,
                     lty = c(rep(NA, length(levels[[1]])), rep(1, length(levels[[2]]))),
                     cex = 0.8, text.font = 2)
            } else if (Faktor %in% fac_names && nf >= 4) {
              error <- "Higher-way interactions cannot be plotted!"
            }
          }
          if (!is.null(error1)) {
            hbox <- RGtk2::gtkHBoxNew()
            vbox$packStart(hbox, FALSE, FALSE, 0)
            label <- RGtk2::gtkLabel(error1)
            hbox$packStart(label, FALSE, FALSE, 0)
          }
          if (!is.null(error)) {
            hbox <- RGtk2::gtkHBoxNew()
            vbox$packStart(hbox, FALSE, FALSE, 0)
            label <- RGtk2::gtkLabel(error)
            hbox$packStart(label, FALSE, FALSE, 0)
          }
        }
        # Create window
        window <- RGtk2::gtkWindow()
        # Add title
        window["title"] <- "Plot"
        # Add a frame
        frame <- RGtk2::gtkFrameNew("Please choose the factor you wish to plot (for interaction type something like group1:group2).")
        window$add(frame)
        # Create vertical container for file name entry
        vbox <- RGtk2::gtkVBoxNew(FALSE, 8)
        vbox$setBorderWidth(24)
        frame$add(vbox)
        # Add horizontal container for every widget line
        hbox <- RGtk2::gtkHBoxNew(FALSE, 8)
        vbox$packStart(hbox, FALSE, FALSE, 0)
        # Add label in first column
        label <- RGtk2::gtkLabelNewWithMnemonic("_Factor")
        hbox$packStart(label, FALSE, FALSE, 0)
        # Add entry in the second column; named "filename"
        filename <- RGtk2::gtkEntryNew()
        filename$setWidthChars(50)
        label$setMnemonicWidget(filename)
        hbox$packStart(filename, FALSE, FALSE, 0)
        # Add an horizontal container to specify parameters
        hbox <- RGtk2::gtkHBoxNew(FALSE, 8)
        vbox$packStart(hbox, FALSE, FALSE, 0)
        label2 <- RGtk2::gtkLabelNewWithMnemonic("_Title")
        hbox$packStart(label2, FALSE, FALSE, 0)
        # Add entry in the second column; named "filename2"
        filename2 <- RGtk2::gtkEntryNew()
        filename2$setWidthChars(10)
        label2$setMnemonicWidget(filename2)
        hbox$packStart(filename2, FALSE, FALSE, 0)
        label3 <- RGtk2::gtkLabelNewWithMnemonic("_lwd")
        hbox$packStart(label3, FALSE, FALSE, 0)
        # Add entry in the second column; named "filename3"
        filename3 <- RGtk2::gtkEntryNew()
        filename3$setWidthChars(10)
        filename3$setText(2)
        label3$setMnemonicWidget(filename3)
        hbox$packStart(filename3, FALSE, FALSE, 0)
        # Add button
        the.buttons <- RGtk2::gtkHButtonBoxNew()
        the.buttons$setBorderWidth(5)
        vbox$add(the.buttons)
        the.buttons$setLayout("spread")
        the.buttons$setSpacing(40)
        buttonOK <- RGtk2::gtkButtonNewFromStock("gtk-ok")
        RGtk2::gSignalConnect(buttonOK, "clicked", plotting)
        the.buttons$packStart(buttonOK, fill=F)
        buttonCancel <- RGtk2::gtkButtonNewFromStock("gtk-close")
        RGtk2::gSignalConnect(buttonCancel, "clicked", window$destroy)
        the.buttons$packStart(buttonCancel, fill=F)
      }
      calculateGUIplot()
    }
    # Output ------------------------------------------------------
    WTS_output <- cbind(WTS_out, WTPS_out)
    colnames(WTS_output) <- cbind ("Test statistic", "df", "p-value",
                                   "p-value WTPS")
    output <- list()
    output$input <- input_list
    output$Descriptive <- descriptive
    output$WTS <- WTS_output
    output$ATS <- ATS_out
  }
  class(output) <- "GFD"
  return(output)
}
