#' A graphical user interface for the package GFD
#' 
#' This function provides a graphical user interface for calculating statistical
#' tests in general factorial designs.
#' 
#' The function produces a GUI for the calculation of the test statistics and 
#' for plotting. Data can be loaded via the "load data" button. The formula, 
#' number of permutations (default: 10,000) and the significance level alpha
#' (default: 0.05) need to be specified. If the plot option is chosen, an
#' additional window opens containing information on the plots.
#' 
#' @export

calculateGUI <- function() {
  ## Run on "Load"
  getDirectory <- function(button, user.data){
    directory <- file.choose()
    RGtk2::gtkEntrySetText(filename, directory)
  }
  ## Run on "OK"
  performStatistics <- function(button, user.data) {
    res <- NULL
    d <- NULL
    error <- NULL
    # Get the information about data and the file
    the.file <- filename$getText()
    the.formula <- formula(filename1$getText())
    the.perm <- as.numeric(filename2$getText())
    the.alpha <- as.numeric(filename3$getText())
    the.plot <- toPlot$active
    the.sep <- sepEntry$getText()
    the.headers <- headersEntry$active
    the.dec <- decEntry$getText()
    d <- read.table(the.file, sep = the.sep, header = the.headers,
                    dec = the.dec)
    res <- GFD(the.formula, d, nperm = the.perm, alpha = the.alpha,
               plot_CI = the.plot, Title = NULL, line_width = 2)
    summary(res)
  }
  # Create window
  window <- RGtk2::gtkWindow()
  # Add title
  window["title"] <- "Tests for general factorial designs"
  # Add a frame
  frame <- RGtk2::gtkFrameNew("Specify data location and formula...")
  window$add(frame)
  # Create vertical container for file name entry
  vbox <- RGtk2::gtkVBoxNew(FALSE, 8)
  vbox$setBorderWidth(24)
  frame$add(vbox)
  # Add horizontal container for every widget line
  hbox <- RGtk2::gtkHBoxNew(FALSE, 8)
  vbox$packStart(hbox, FALSE, FALSE, 0)
  # Add label in first column
  label <- RGtk2::gtkLabelNewWithMnemonic("_File name")
  hbox$packStart(label, FALSE, FALSE, 0)
  # Add entry in the second column; named "filename"
  filename <- RGtk2::gtkEntryNew()
  filename$setWidthChars(50)
  label$setMnemonicWidget(filename)
  hbox$packStart(filename, FALSE, FALSE, 0)
  # Add label in first column
  label1 <- RGtk2::gtkLabelNewWithMnemonic("_Formula")
  hbox$packStart(label1, FALSE, FALSE, 0)
  # Add entry in the second column; named "filename1"
  filename1 <- RGtk2::gtkEntryNew()
  filename1$setWidthChars(50)
  label1$setMnemonicWidget(filename1)
  hbox$packStart(filename1, FALSE, FALSE, 0)
  # Add an horizontal container to specify parameters
  hbox <- RGtk2::gtkHBoxNew(FALSE, 8)
  vbox$packStart(hbox, FALSE, FALSE, 0)
  label2 <- RGtk2::gtkLabelNewWithMnemonic("_nperm")
  hbox$packStart(label2, FALSE, FALSE, 0)
  # Add entry in the second column; named "filename2"
  filename2 <- RGtk2::gtkEntryNew()
  filename2$setWidthChars(10)
  filename2$setText(10000)
  label2$setMnemonicWidget(filename2)
  hbox$packStart(filename2, FALSE, FALSE, 0)
  label3 <- RGtk2::gtkLabelNewWithMnemonic("_alpha")
  hbox$packStart(label3, FALSE, FALSE, 0)
  # Add entry in the second column; named "filename3"
  filename3 <- RGtk2::gtkEntryNew()
  filename3$setWidthChars(10)
  filename3$setText(0.05)
  label3$setMnemonicWidget(filename3)
  hbox$packStart(filename3, FALSE, FALSE, 0)
  # Add an horizontal container to specify input file options
  # are headers included in the file?
  hbox <- RGtk2::gtkHBoxNew(FALSE, 8)
  vbox$packStart(hbox, FALSE, FALSE, 0)
  label <- RGtk2::gtkLabelNewWithMnemonic("_Headers?")
  hbox$packStart(label, FALSE, FALSE, 0)
  headersEntry <- RGtk2::gtkCheckButton()
  headersEntry$active <- TRUE
  hbox$packStart(headersEntry, FALSE, FALSE, 0)
  label$setMnemonicWidget(headersEntry)
  # what separator is used?
  label <- RGtk2::gtkLabelNewWithMnemonic("Col. _Separator?")
  hbox$packStart(label, FALSE, FALSE, 0)
  sepEntry <- RGtk2::gtkEntryNew()
  sepEntry$setWidthChars(1)
  sepEntry$setText("")
  hbox$packStart(sepEntry, FALSE, FALSE, 0)
  label$setMnemonicWidget(sepEntry)
  # what's the character used for decimal points?
  label <- RGtk2::gtkLabelNewWithMnemonic("_Dec. character?")
  hbox$packStart(label, FALSE, FALSE, 0)
  decEntry <- RGtk2::gtkEntryNew()
  decEntry$setWidthChars(1)
  decEntry$setText(".")
  hbox$packStart(decEntry, FALSE, FALSE, 0)
  label$setMnemonicWidget(decEntry)
  # Add separator
  vbox$packStart(RGtk2::gtkHSeparatorNew(), FALSE, FALSE, 0)
  # Add plot-option
  hbox <- RGtk2::gtkHBoxNew(FALSE, 8)
  vbox$packStart(hbox, FALSE, FALSE, 0)
  label <- RGtk2::gtkLabelNewWithMnemonic("Plot _Results?")
  hbox$packStart(label, FALSE, FALSE, 0)
  toPlot <- RGtk2::gtkCheckButton()
  hbox$packStart(toPlot, FALSE, FALSE, 0)
  # Add button
  the.buttons <- RGtk2::gtkHButtonBoxNew()
  the.buttons$setBorderWidth(5)
  vbox$add(the.buttons)
  the.buttons$setLayout("spread")
  the.buttons$setSpacing(40)
  buttonOK <- RGtk2::gtkButtonNewFromStock("gtk-ok")
  buttonLoad <- RGtk2::gtkButtonNewFromStock("load data")
  RGtk2::gSignalConnect(buttonOK, "clicked", performStatistics)
  RGtk2::gSignalConnect(buttonLoad, "clicked", getDirectory)
  the.buttons$packStart(buttonOK, fill=F)
  the.buttons$packStart(buttonLoad, fill=F)
  buttonCancel <- RGtk2::gtkButtonNewFromStock("gtk-close")
  RGtk2::gSignalConnect(buttonCancel, "clicked", window$destroy)
  the.buttons$packStart(buttonCancel, fill=F)
}
