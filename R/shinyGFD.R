#' A shiny app for the package GFD
#'
#' This function provides a shiny app for calculating
#' GFD and QANOVA test statistics and respective p-values.
#'
#'@author Philipp Steinhauer
#'
#'
#' @aliases GFDsurvGUI
#'
#'
#' @import shiny
#' @import utils
#' @importFrom tippy tippy_this
#' @importFrom shinythemes shinytheme
#' @importFrom shinyWidgets numericInputIcon
#' @importFrom shinyjs useShinyjs hide hidden show
#' @importFrom stats as.formula qnorm pbinom bw.nrd0
#' @export

GFD_GUI <- function() {
  requireNamespace("shiny", quietly = TRUE)
  
  if (!("package:shiny" %in% search())) {
    attachNamespace("shiny")
  }
  requireNamespace("shinyWidgets", quietly = TRUE)
  if (!("package:shinyWidgets" %in% search())) {
    attachNamespace("shinyWidgets")
  }
  requireNamespace("tippy", quietly = TRUE)
  
  if (!("package:tippy" %in% search())) {
    attachNamespace("tippy")
  }
  
  
  ui <- fluidPage(theme = shinythemes::shinytheme("cerulean"),
                  shinyjs::useShinyjs(),
                  titlePanel("Tests for General Factorial Designs"),
                  sidebarLayout(
                    sidebarPanel(
                      splitLayout(cellArgs = list(style = "vertical-align: middle"),
                        fileInput("infile", "Choose CSV File",
                                  accept = c(
                                    "text/csv",
                                    "text/comma-separated-values,text/plain",
                                    ".csv")),
                        checkboxInput("header", "Header", TRUE),
                        selectInput("sep","Seperator in csv", c(",",
                                                                ";",
                                                                ".",
                                                                "|"))
                      ),
                      
                      tags$head(tags$style(HTML("
                              .shiny-split-layout > div {
                                overflow: visible;
                              }
                              "))), #for selectinput in splitlayout with full dropdown view
                      tags$style(HTML("
                                 input[type=number] {
                                                              -moz-appearance:textfield;
                                                    }
                                  input[type=number]::{
                                                  -moz-appearance:textfield;
                                                    }
                        input[type=number]::-webkit-outer-spin-button,
                        input[type=number]::-webkit-inner-spin-button {
                        -webkit-appearance: none;
                        margin: 0;
                        }
                        ")),
                      
                      
                      h3(id="titleLoadData","Load dataset first!", style = "color:red"),
                      
                      shinyjs::hidden(
                        selectInput("Method", "Select Testing Method:",
                                    c("GFD: Mean-based tests for General Factorial Designs" = "GFD",
                                      "QANOVA: Quantile-based tests for General Factorial Designs"= "QANOVA"))
                      ),
                      
                      splitLayout(cellWidths = c("50%","5%","45%"),cellArgs = list(style = "vertical-align: sub"),
                                  shinyjs::hidden(
                                    textInput("formula", "Formula ", "Response  ~ FactorA * FactorB")
                                  ),
                                  shinyjs::hidden(
                                    actionButton("infoButton", "", icon = icon("info-circle"))
                                  ),
                                  
                                  tippy::tippy_this("infoButton", "Interaction effects need to be specified! Example:<br><br>
                                              - 1 Factor: <br>
                                              response  ~ factorA <br><br>
                                              - 2 Factors without interactions:
                                              <br> response  ~ factorA + factorB
                                              <br>  <br>- 2 Factors with interactions:
                                              <br> response  ~ factorA + factorB + factorA:factorB
                                              <br> or
                                              <br> response  ~ factorA * factorB "
                                                    , placement = "right")
                      ),
                      
                      splitLayout(cellWidths = c("60%"),
                                  shinyjs::hidden(
                                    checkboxInput("nested", "Are the levels of nested factors the same for each level main factor?", FALSE)
                                  )
                      ),
                      
                      splitLayout(cellWidths = c("35%","60%","5%"),
                                  # style = "border: 1px solid silver;",
                                  # cellArgs = list(style = "padding: 6px"),
                                  shinyjs::hidden(
                                  checkboxGroupInput("median", "Single quantiles",selected = c("median"),
                                                     choiceNames = list("Median"),
                                                     choiceValues = list("median"))
                                  ),
                                  shinyjs::hidden(
                                    selectInput("singles","user-specified quantiles",
                                                paste0(seq(5,95,5),"%"),
                                                multiple=TRUE,selectize = TRUE)
                                  )
                      ),
                      
                      
                    
                      splitLayout(cellWidths = c("35%","13%","5%","13%"),
                                   # style = "border: 1px solid silver;", 
                                   # cellArgs = list(style = "padding: 6px"),
                                  shinyjs::hidden(
                                    checkboxGroupInput("ranges", "Range between two quantiles",
                                                       choiceNames = list("IQR (75% - 25%)","90% - 10%"),
                                                       choiceValues = list("IQR","diff_90_10"))
                                  ),
                                  shinyjs::hidden(
                                    shinyWidgets::numericInputIcon(
                                    inputId = "range_user",
                                    label = "user-specified ranges",
                                    value = 95,
                                    min = 0,
                                     max = 100,
                                    icon = list(NULL, icon("percent"))
                                    )
                                  ),
                                  
                                  
                                  
                                  shinyjs::hidden(
                                                 h2(id="title_range","  -", style = "color:black,font-size:500px;"),
                                                 tags$style(type='text/css', "#title_range {margin-top: 25px;}")
                                  ),
                                    
                                  shinyjs::hidden(
                                    shinyWidgets::numericInputIcon(
                                      inputId = "range_user2",
                                      label = ".",
                                      value = 5,
                                      min = 0,
                                      max = 100,
                                      icon = list(NULL, icon("percent"))
                                    )
                                    #tags$style(type='text/css', "#range_user2 {margin-top: 5px;}")
                                  )
                      ),   
                      
                      
                            
                      
                      shinyjs::hidden(
                        checkboxInput("plots", "Plot", FALSE)
                      ),
                      
                      splitLayout(cellWidths = c("60%","5%"),cellArgs = list(style = "vertical-align: sub"),
                                  
                                  shinyjs::hidden(
                                    textInput("plot_factor", "Please choose the factor you wish to plot","Name of factor")
                                  ),
                                  shinyjs::hidden(
                                    actionButton("infoButton2", "", icon = icon("info-circle"))
                                  ),
                                  
                                  tippy::tippy_this("infoButton2", "Interaction effects need to be specified! Example:<br><br>
                                              - 1 Factor: <br>
                                              factorA <br><br>
                                              - 2 Factors with interactions:
                                              <br> factorA:factorB"
                                                    , placement = "right")
                      ),
                      
                      
                      
                      splitLayout(cellWidths = c("50%","10%","0%","30%"),
                                  shinyjs::hidden(
                                    selectInput("Var_method", "Confidence intervals based on",
                                              c("t-quantiles" = "t-quantile",
                                                "permutation distribution" = "perm"))
                                  ),
                                  shinyjs::hidden(
                                    actionButton("infoButton3", "", icon = icon("info-circle"))
                                  ),
                                  
                                  tippy::tippy_this("infoButton3", "Method for the variance estimation of the sample quantiles
                                  can be chosen:<br><br>
                                              - the interval-based estimator of Price and Bonett (2001) <br><br>
                                              - Efron's bootstrap method  <br><br> - kernel density approach"
                                                    , placement = "right"),
                                  shinyjs::hidden(
                                    numericInput("alpha2", "confidence level", value = 0.95,min = 0, max = 1)
                                  ),
                                  tags$style(type='text/css', "#infoButton3 {margin-top: 25px;}")
                      ),
                      
                     
                      splitLayout(cellWidths = c("30%","15%"),
                                  shinyjs::hidden(
                                    numericInput("nperm", "Number of permutations", value = 1999)
                                  )
                      ),
                      
                     
                      
                      shinyjs::hidden(
                        actionButton("process", "Calculate", class = "btn-primary")
                      )
                      
                     
                      , width = 6
                    ),
                    
                  
                    
                    mainPanel(
                      
                      verbatimTextOutput("result"),
                      plotOutput("result_plot"),
                      width = 6
                      
                    )
                  )
  )
  
  
  server <- function(input, output,session) {
    
    datasetInput <- reactive({
      
      req(input$infile)
      
      if (is.null(input$infile))
        return(NULL)
      read.csv(input$infile$datapath, header = input$header, sep = as.character(input$sep))
    })
    
    
    observeEvent(input$infile, {
      
      if(is.null(input$infile)){
        shinyjs::hide(id = "Method")
        shinyjs::hide(id = "nperm")
        shinyjs::hide(id = "Var_method")
        shinyjs::hide(id = "formula")
        shinyjs::hide(id = "infoButton")
        shinyjs::hide(id = "alpha2")
        shinyjs::hide(id = "process")
        shinyjs::hide(id = "singles")
        shinyjs::hide(id = "median")
        shinyjs::hide(id = "ranges")
        shinyjs::hide(id = "range_user")
        shinyjs::hide(id = "range_user2")
        shinyjs::hide(id = "plots")
        shinyjs::hide(id = "plot_factor")
        shinyjs::hide(id = "infoButton2")
        shinyjs::hide(id = "infoButton3")
        shinyjs::hide(id = "title_range")
        shinyjs::hide(id = "title_quantiles")
        
      }else{
        shinyjs::show(id = "Method")
        shinyjs::show(id = "formula")
        shinyjs::show(id = "infoButton")
        shinyjs::show(id = "nperm")
        shinyjs::show(id = "Var_method")
        shinyjs::show(id = "nested")
        shinyjs::hide(id = "titleLoadData")
        shinyjs::show(id = "process")
        shinyjs::show(id = "alpha2")
        
        
        observeEvent(input$Method, {
          
          if (input$Method == "GFD") {
            
            shinyjs::hide(id = "singles")
            shinyjs::hide(id = "median")
            shinyjs::hide(id = "ranges")
            shinyjs::hide(id = "range_user")
            shinyjs::hide(id = "range_user2")
            shinyjs::hide(id = "title_quantiles")
            shinyjs::hide(id = "title_range")
            shinyjs::show(id = "plots")
            shinyjs::hide(id = "infoButton3")
            
            
            observeEvent(input$plots, {
              if (input$plots) {
                
                shinyjs::show(id = "plot_factor")
                shinyjs::show(id = "infoButton2")
                updateSelectInput(session, "Var_method","Confidence intervals based on",
                                  c("t-quantiles" = "t-quantile",
                                    "permutation distribution" = "perm"))
                shinyjs::show(id = "Var_method")
                shinyjs::show(id = "alpha2")
                updateNumericInput(session,"alpha2", "confidence level", value = 0.95,min = 0, max = 1)
                shinyjs::show(id = "alpha2")
                
              } else{
                
                shinyjs::hide(id = "plot_factor")
                shinyjs::hide(id = "infoButton2")
                shinyjs::hide(id = "Var_method")
                shinyjs::hide(id = "alpha2")
                
              }
              
            })
            
            
          }
          
          
          if (input$Method == "QANOVA") {
            
            updateSelectInput(session, "Var_method","Variance estimation of the sample quantiles by:",
                              c("Interval-based estimator" = "interval",
                                "Bootstrap method" = "boot",
                                "Kernel density approach " = "kernel"))
            shinyjs::show(id = "singles")
            shinyjs::show(id = "median")
            shinyjs::show(id = "ranges")
            shinyjs::show(id = "range_user")
            shinyjs::show(id = "range_user2")
            shinyjs::hide(id = "plots")
            shinyjs::hide(id = "plot_factor")
            shinyjs::hide(id = "infoButton2")
            shinyjs::show(id = "infoButton3")
            shinyjs::show(id = "Var_method")
            shinyjs::show(id = "title_quantiles")
            shinyjs::show(id = "title_range")
            
            observeEvent(input$Var_method, {
              if (input$Var_method == "interval") {
                updateNumericInput(session,"alpha2", "Confidence level of Interval", value = 0.95,min = 0, max = 1)
                shinyjs::show(id = "alpha2")
                
              }
              
              if (input$Var_method == "boot" || input$Var_method == "kernel"){ 
                shinyjs::hide(id = "alpha2")
              }
            
              if (input$Var_method == "t-quantile" || input$Var_method == "perm") {
                updateNumericInput(session,"alpha2", "confidence level", value = 0.95,min = 0, max = 1)
                shinyjs::show(id = "alpha2")
                
              }
            
              
            })
          }
        
          
        })## observeevent
        
        
        
      }
    })
    
    
    
    
    
    
    observeEvent(input$process, {
      
      if (input$formula == "Response  ~ FactorA*FactorB") {
        
        output$result <- renderPrint({
          "'formula' missing or invalid"
        })
        
      } else {
          
          if (input$Method == "GFD" ){
            data <- as.data.frame(datasetInput())
            
            showModal(modalDialog("Calculating!"))
            output_GFD <- GFD(formula= as.formula(isolate(input$formula)),
                                   data = isolate(data),
                                   nperm = isolate(input$nperm),
                                   alpha =  1-isolate(input$alpha2),
                                   nested.levels.unique = isolate(input$nested),
                                   CI.method = isolate(input$Var_method)
                              )
            removeModal()
            
            output$result <- renderPrint({
              output_GFD
            })
            
            
             if(input$plots){
               #Plotting Function aus utility uebernommen
               object <- output_GFD
               
               a <- object$plotting
               b <- object$Descriptive
               fac.names <- a$fac_names
               exist <- hasArg(factor) 
               
              Faktor <- input$plot_factor
               x.label <- input$plot_factor
               
               match.arg(Faktor, fac.names)
               
               # default values
               args <- list(plot.object = a, descr.object = b, factor = Faktor,
                            lwd =2, ylab = "Means", xlab = x.label, col = 1:length(fac.names), pch = 1:18, legendpos = "topright")
               

               
               output$result_plot <- renderPlot({
                 do.call(plotting, args = args)
                 
               })
               
             }
          }
        
        if (input$Method == "QANOVA" ){
          data <- as.data.frame(datasetInput())
          
         
          givenQuant <- unlist(list(isolate(input$median),isolate(input$ranges)))
          extraQuant <- unlist(list(isolate(input$singles)))
          range <- unlist(list(isolate(input$range_user),isolate(input$range_user2)))
          all_null <- unlist(list(givenQuant,extraQuant,range))
          # output$result <- renderPrint({
          #          list(givenQuant,extraQuant,range,all_null, sum(is.na(range))>0)
          # 
          #     })

          
          if(is.null(givenQuant) && is.null(extraQuant) && sum(is.na(range))>0){

          output$result <- renderPrint({
            "No quantiles had been choosen"

          })} else{
            showModal(modalDialog("Calculating!"))


            if(is.null(givenQuant)){
              matrix <- matrix(9999)
              quantiles <- c(9999)
            } else{
              quant <- c("median", "IQR","diff_90_10")
              quant_TF <- quant %in% givenQuant
              matrix <- matrix(c(1,0,0,0,0,0,1,-1,0,0,0,0,0,1,-1),3,byrow =T)
              matrix <- matrix[quant_TF,]
              quantiles <- c(0.5,0.75,0.25,0.9,0.1)
                if(sum(quant_TF)>1){
                  quantiles <- quantiles[which(apply(abs(matrix),2,sum)>0)]

                  matrix <- matrix[,which(apply(abs(matrix),2,sum)>0)]
                  } else {
                  quantiles <- quantiles[which(abs(matrix)==1)]
                  matrix  <-  t(as.matrix(matrix[which(abs(matrix)==1)]))
                }
            }

            if(!is.null(extraQuant)){

            quantiles_single <- as.numeric(sub("%", "", extraQuant))/100
            quan_single <- quantiles_single %in% quantiles
            for(i in 1:length(quan_single)){
              if(quan_single[i]){
                vector <- rep(0,ncol(matrix))
                vector[which(quantiles==quantiles_single[i])] <- 1
                matrix <- rbind(matrix,vector)
              }else{
                matrix <- cbind(matrix,0)
                matrix <- rbind(matrix,c(rep(0,ncol(matrix)-1),1))
              }
            }
            quantiles <- unique(c(quantiles,quantiles_single))
            }


            if(sum(is.na(range))==0){
              range_single <- range/100
              matrix <- rbind(matrix,0)
              if(range_single[1] %in% quantiles){
                matrix[nrow(matrix),which(quantiles==range_single[1])] <- 1
              }else{
                matrix <- cbind(matrix,0)
                matrix[nrow(matrix),ncol(matrix)] <- 1
              }

              if(range_single[2] %in% quantiles){
                matrix[nrow(matrix),which(quantiles==range_single[2])] <- -1
              }else{
                matrix <- cbind(matrix,0)
                matrix[nrow(matrix),ncol(matrix)] <- -1
              }
              quantiles <- unique(c(quantiles,range_single))

            }
            if(9999 %in% quantiles){
              quantiles <- quantiles[-1]
              matrix <- matrix[-1,-1]
            }
            if(!is.null(nrow(matrix))){
            rownames(matrix) <- 1:nrow(matrix)
            matrix <- matrix[,order(quantiles)]
            } else{
              matrix <- matrix(matrix,1)
              matrix <- matrix(matrix[,order(quantiles)],1)
            }
            quantiles <- sort(quantiles)
            
            output_qanova <- QANOVA(formula= as.formula(isolate(input$formula)),
                               data = isolate(data),
                               quantiles = quantiles,
                               lin_mat = matrix,
                               var_method = isolate(input$Var_method),
                               nperm = isolate(input$nperm),
                               var_level = isolate(input$alpha2),
                               nested.levels.unique = isolate(input$nested)
             )
            removeModal()

            output$result <- renderPrint({
                output_qanova

           })

          }
        }
         
          
       }
        
    }) #end of observeEvent(input$process
    
  }
  
  
  shinyApp(ui = ui, server = server)
  
}

