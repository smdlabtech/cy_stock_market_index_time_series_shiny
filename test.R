library(shiny)
library(PerformanceAnalytics)
library(quantmod)
library(rugarch)


#----------------------------------#
# UI : user Interface
#----------------------------------#

ui <- fluidPage(
  titlePanel("Analyse des rendements financiers"),
  sidebarLayout(
    sidebarPanel(
      selectInput("symbol1", "Symbole 1:", choices = c("AMZN", "FB")),
      selectInput("symbol2", "Symbole 2:", choices = c("AMZN", "FB")),
      dateRangeInput("dateRange", "Période :", start = "2012-06-01", end = "2021-12-01"),
      actionButton("analyze", "Analyser")
    ),
    mainPanel(
      plotOutput("priceplot", height = "300px"),
      plotOutput("returnplot", height = "300px"),
      plotOutput("volatilityplot", height = "300px"),
      plotOutput("acfplot", height = "600px"),
      tableOutput("statsTable"),
      plotOutput("infocriteriaplot", height = "400px")
    )
  )
)


#----------------------------------#
# SERVER
#----------------------------------#
server <- function(input, output) {
  
  # Charger les données
  data <- reactive({
    getSymbols(c(input$symbol1, input$symbol2), from = input$dateRange[1], to = input$dateRange[2])
    data <- cbind(get(input$symbol1)[, "Adjusted"], get(input$symbol2)[, "Adjusted"])
    colnames(data) <- c(input$symbol1, input$symbol2)
    data
  })
  
  # Calculer les rendements logarithmiques
  ret <- reactive({
    ret1 <- CalculateReturns(data()[, 1], method = "log")
    ret2 <- CalculateReturns(data()[, 2], method = "log")
    cbind(ret1, ret2)
  })
  
  # Afficher les graphiques
  output$priceplot <- renderPlot({
    plot.zoo(data(), main = "Évolution des cours")
  })
  
  output$returnplot <- renderPlot({
    plot.zoo(ret(), main = "Évolution des rendements")
  })
  
  output$volatilityplot <- renderPlot({
    dataToPlot <- cbind(ret()[, 1]^2, ret()[, 2]^2)
    colnames(dataToPlot) <- c(paste0("abs(", input$symbol1, ".ret)"), paste0("abs(", input$symbol2, ".ret)"))
    plot.zoo(dataToPlot, main = "Évolution de la volatilité")
  })
  
  output$acfplot <- renderPlot({
    par(mfrow = c(3, 2))
    acf(ret()[, 1], main = paste0(input$symbol1, " Rendements"))
    acf(ret()[, 2], main = paste0(input$symbol1, " Rendements"))
    acf(ret()[, 1]^2, main = paste0(input$symbol1, " Rendements^2"))
    acf(ret()[, 2]^2, main = paste0(input$symbol2, " Rendements^2"))
    acf(abs(ret()[, 1]), main = paste0(input$symbol1, " abs(Rendements)"))
    acf(abs(ret()[, 2]), main = paste0(input$symbol2, " abs(Rendements)"))
  })
  
  output$statsTable <- renderTable({
    A <- table.Stats(ret()[, 1])
    B <- table.Stats(ret()[, 2])
    tab_desc <- cbind(A, B)
    tab_desc
  })
  
  output$infocriteriaplot <- renderPlot({
    arch11.spec <- ugarchspec(variance.model = list(garchOrder = c(1, 1)), mean.model = list(armaOrder = c(0, 0)))
    AMZN.arch11.fit <- ugarchfit(spec = arch11.spec, data = ret()[, 1])
    FB.arch11.fit <- ugarchfit(spec = arch11.spec, data = ret()[, 2])
    
    model.list <- list(garch11 = AMZN.arch11.fit, egarch11 = AMZN.egarch11.fit, gjrgarch11 = AMZN.gjrgarch11.fit, aparch11.1 = AMZN.aparch11.1.fit)
    info.mat_amzn <- sapply(model.list, infocriteria)
    rownames(info.mat_amzn) <- rownames(infocriteria(AMZN.arch11.fit))
    info.mat_amzn <- t(info.mat_amzn)
    info.mat_amzn <- abs(info.mat_amzn)
    
    x <- c(1, 2, 3, 4)
    plot((info.mat_amzn[1, 1:4]) ~ x, xlab = "Critères", ylab = "Valeur des critères", type = 'l', col = "blue", main = "Critères d'information", xlim = c(1, 4), ylim = c(5.195, 5.25))
    lines(info.mat_amzn[2, 1:4], col = "gray")
    lines(info.mat_amzn[3, 1:4], col = "black")
    lines(info.mat_amzn[4, 1:4], col = "orange")
    legend("topleft", legend = c("garch(1,1)", "egarch(1,1)", "gjrgarch(1,1)", "aparch(1,1)"), fill = c("blue", "gray", "black", "orange"))
  })
  
  observeEvent(input$analyze, {
    # Effectuer l'analyse
  })
}


#----------------------------------#
# LAUCH APP
#----------------------------------#
shinyApp(ui, server)

