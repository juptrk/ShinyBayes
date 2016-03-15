################## IMPORTS, FIXED VARIABLES ##########

library(BEST)
library(DT)
library(prob)
library(BayesFactor)
library(MASS)
library(ggplot2)
library(rjags)
library(ggmcmc)
library(dplyr)
library(stats)

aComp = 1
bComp = 1
alpha = 0.05
threshold = 0.001


# Function for getting the meaning of the bf
get_table_content = function(bf){
  if (bf <= 1) {
    result = data.frame(c("1"), c("irrelevant data"))
  }
  else if (bf <= 3) {
    result = data.frame(c("1 - 3"), c("hardly worth ink or breath"))
  }
  else if (bf <= 6) {
    result = data.frame(c("3 - 6"), c("anecdotal"))
  }
  else if (bf <= 10) {
    result = data.frame(c("6 - 10"), c("now we're talking: substantial"))
  }
  else if (bf <= 30) {
    result = data.frame(c("10 - 30"), c("strong"))
  }
  else if (bf <= 100) {
    result = data.frame(c("30 - 100"), c("very strong"))
  }
  else {
    result = data.frame(c("100+"), c("decisive"))   
  }
}

# Function for fix N
BinomForFixN = function(k, N, theta){
  result = choose(N, k) *  (theta^k) * (1- theta)^(N-k)
}

# Function for fix k
BinomForFixK = function(N, k, theta){
  if (N != 0){
    result = (k/N) * choose(N, k) *  (theta^k) * (1- theta)^(N-k)
  }
  else {
    result = 0
  }
}

# Function for plot beta
pl.beta <- function(a,b,l,u,theta, asp = if(isLim) 1, ylim = if(isLim) c(0,1.1)){
  
  if(isLim <- a == 0 || b == 0 || a == Inf || b == Inf){
    eps <- 1e-10
    x <- c(0, eps, (1:7)/16, 1/2+c(-eps,0,eps), (9:15)/16, 1-eps,1)
  } 
  else{
    x <- seq(0, 1, length = 1025)
  }
  
  p <- ggplot() + geom_line(aes(x,dbeta(x,a,b))) 
  
  p <- p + ylab("beta(a,b)") + ggtitle(paste("beta(",a,",",b,")"))
  
  p <- p + geom_segment(aes(x = l, y = dbeta(l,a,b),
                            xend = u, yend = dbeta(u,a,b),
                            colour = "hdi"))
  
  p <- p + geom_point(aes(x = l, y = dbeta(l,a,b)), col="red")
  
  p <- p + geom_text(aes(x=l-0.05, y=dbeta(l,a,b)+.12,
                         label=paste(round(l, 3))), col="red")
  
  p <- p + geom_text(aes(x=u+0.05, y=dbeta(u,a,b)+.12,
                         label=paste(round(u, 3))), col="red")
  
  p <- p + geom_point(aes(x = u, y = dbeta(u,a,b)), col="red")
  
  p <- p + geom_segment(aes(x=theta, y = 0, xend=theta, yend = 3), col="blue")
  
  p <- p + geom_point(aes(x=theta, y = 3), col="blue")
  
  p <- p + geom_text(aes(x=theta+0.05,y=3,label=paste(theta)), col="blue")
  
  p
}

#From Kruschke p.270
pD = function(k,N,a,b) { beta(k+a,N-k+b) / beta(a,b) }


# Define server logic required to draw a histogram
shinyServer(function(input, output, session) {
  observe({
    updateSliderInput(
      session,
      "k",  
      value = input$k,
      min = ifelse((input$option == "k fix"), 1, 0),
      max = input$N)
  })
  
  # Calculates hdi
  hdi_coins <- reactive({
    hdi_size = input$hdi / 100
    
    # posterior
    betaParam1 = aComp + input$k
    betaParam2 = bComp + input$N - input$k
    
    # hdi
    HDI = hdi(qbeta, hdi_size, shape1 = betaParam1, shape2 = betaParam2)
    
    list(lower = HDI[1],
         upper = HDI[2],
         beta1 = betaParam1,
         beta2 = betaParam2)
  })
  
  # calculates the p value
  pvalue = reactive({
    k = input$k
    N = input$N
    theta = input$theta
    values <- NULL
    pvalue <- 0
    plotValues <- NULL
    
    if (input$option == "N fix") {
      
      valueFork <- BinomForFixN(k, N, theta)
      max <- valueFork
      
      for(i in 1:N+1){
        values[i] <- BinomForFixN(i-1, N, theta)
        if(values[i] > max){
          max = values[i]
        }
        if(values[i] <= valueFork){
          pvalue <- pvalue + values[i]
        }
      }
      
      values[1] = values[N+1]
      
      xachse = seq(from=0, to=N)
      
      list(values = values,
           valueFork = valueFork,
           pvalue = pvalue,
           max = max,
           xachse = xachse)
    }
    else if (input$option == "k fix") {
      valueForN <- BinomForFixK(N, k, theta)
      max <- valueForN
      values[1] = BinomForFixK(k,k, theta)
      plotValues[1] = values[1]
      
      if(values[1] <= valueForN){
        pvalue = pvalue + values[1]
      }
      if(values[1] > max) {
        max = values[1]
      }
      
      index <- k
      if (k == 0) {
        index = 1
      }
      
      for(i in (k+1):(index*300)){
        values[i-k+1] <- BinomForFixK(i, k, theta)
        if(values[i-k+1] > threshold){
          plotValues[i-k+1] = values[i-k+1]
        }
        else if(values[i-k+1] >= threshold && i < (N+10)){
          plotValues[i-k+1] = values[i-k+1]
        }
        
        if(values[i-k+1] > max){
          max = values[i-k+1]
        }
        
        if(values[i-k+1] <= valueForN){
          pvalue <- pvalue + values[i-k+1]
        }
      }
      
      xachse = seq(from=k, to= k+ length(plotValues)-1)
      list(valueForN = valueForN,
           pvalue = pvalue,
           max = max,
           xachse = xachse,
           values = plotValues)
    }
  })
  
  bf_coins <- reactive({
    theta = input$theta
    k = input$k
    N = input$N
    BF = (theta^k * (1- theta)^(N-k)) / pD(k,N,aComp,bComp)
    list(BF = BF)
  })
  
  corr_set = reactive({
    
    n <- input$N_cor
    corr <- input$corr
    
    data = mvrnorm(n,
                   mu = c(0,0),
                   Sigma = matrix(c(1,corr,corr,1), ncol = 2),
                   empirical = TRUE)
    
    colnames(data) <- c("x", "y")
    
    data.frame(data)
  })
  
  linear_bf = reactive({
    
    data = corr_set()
    
    bf <- lmBF(y ~ x, data = data, iterations = 10000)
    
  })
  
  output$plot_estimation <- renderPlot({
    
    beta1 = hdi_coins()$beta1
    beta2 = hdi_coins()$beta2
    lower = hdi_coins()$lower
    upper = hdi_coins()$upper
    
    pl.beta(beta1, beta2, lower, upper, input$theta) 
  })
  
  
  output$text_rope <- renderText({
    
    lower = hdi_coins()$lower
    upper = hdi_coins()$upper
    
    lowerBoundaryRope = input$theta - input$ROPE
    upperBoundaryRope = input$theta + input$ROPE
    boolean = ifelse((lowerBoundaryRope >= lower && upperBoundaryRope <= upper), TRUE, FALSE)
    
    if (boolean) {
      paste("ROPE with epsilon = ", input$ROPE, " confirmes the HDI result.")
    }
    else {
      paste("ROPE with epsilon = ", input$ROPE, " does not confirm the HDI result.")
    }
    
  })
  
  output$plot_pvalue <- renderPlot({
    k = input$k
    N = input$N
    theta = input$theta
    
    if (input$option == "N fix") {
      values = pvalue()$values
      valueFork = pvalue()$valueFork
      pvalue = pvalue()$pvalue
      max = pvalue()$max
      xachse = pvalue()$xachse
      
      barplot (height = values,col = ifelse(values <= valueFork, 'deeppink', 'green2'),
               xlim = c(0,N+2), xlab = "k", names.arg = xachse,
               ylim = c(0, max+0.05),space = 0, axes = TRUE,
               axisnames = TRUE, axis.lty = 0.2,
               ylab=paste(paste("P( k | N = " , N, ""),
                          paste(paste(",", theta, ""), ")", ""), ""))
      segments(0,valueFork, N+1, valueFork, lwd=3, col='red')
      text(N-2, valueFork+.01, paste("p = ", round(pvalue, 4)),
           col = "red", cex = 1 )
    }
    
    else if (input$option == "k fix") {
      valueForN = pvalue()$valueForN
      pvalue = pvalue()$pvalue
      max = pvalue()$max
      xachse = pvalue()$xachse
      plotValues = pvalue()$values
      
      barplot (height = plotValues,col = ifelse(plotValues <= valueForN,
                                                'deeppink', 'green2'),
               xlim = c(0, length(plotValues)-1), xlab = "N",
               names.arg = xachse, ylim = c(0, max+0.05),space = 0,
               axes = TRUE, axisnames = TRUE, axis.lty = 0.2,
               ylab=paste(paste("P( N | k = " , k, ""),
                          paste(paste(",", theta, ""), ")", ""), ""))
      segments(0,valueForN,length(plotValues), valueForN, lwd=3, col='red')
      text(length(plotValues) - 2, valueForN+.01, paste("p = ", round(pvalue, 4)),
           col = "red", cex = 1 )
      
    }
    else {
      paste("Sorry, something went wrong with the slider input!")
    }
  })
  
  output$text_pvalue <- renderText({
    k = input$k
    N = input$N
    theta = input$theta
    values <- NULL
    pvalue <- 0
    
    #### TODO: k = 1 und k = 0 abfangen schmeißt fehler!
    if (input$option == "N fix") {
      pvalue <- pvalue()$pvalue
      if (pvalue <= alpha) {
        paste("The p-value is ", round(pvalue, 4),
              ". Thus it is significant with alpha = .05.")
      }
      else {
        paste("The p-value is ", round(pvalue, 4),
              ". Thus it is not significant with alpha = .05.")
      }
    }
    else if (input$option == "k fix") {
      pvalue <- pvalue()$pvalue
      
      if (pvalue <= alpha) {
        paste("The p-value is ", round(pvalue, 4),
              ". Thus it is significant with alpha = .05.")
      }
      else {
        paste("The p-value is ", round(pvalue, 4),
              ". Thus it is not significant with alpha = .05.")
      }
    }
  })
  

  output$text_comparison <- renderText({
    BF = bf_coins()$BF
    paste("The Bayes factor is ", round(BF, 4), ".")
  })
  
  output$table_bf <- renderDataTable({
    BF = bf_coins()$BF
    
    bf_frame = get_table_content(BF)
    
    colnames(bf_frame) <- c("BF (M0 > M1)", "interpretation")
    bf_frame
  }, rownames = FALSE, options = list(paging = FALSE, searching = FALSE))
  
  
  output$cor_plot <- renderPlot({
    
    data = corr_set()
    
    p <- ggplot(data = data, aes(x=data$x, y=data$y)) + xlim(-3,3) + ylim(-3,3)
    
    p <- p + geom_point(data= data,aes(x=data$x, y=data$y),color='red')
    
    p <- p + xlab(input$Xaxis) + ylab(input$Yaxis) 
    
    p <- p + ggtitle("Plot for the selected correlation")
    
    p <- p + geom_smooth(method = "lm", colour = "black", level = 0)
    
    p
  })
  
  
  output$linear_bayes_text <- renderText({
    bf <- extractBF(linear_bf())
    
    paste("The Bayes factor is ", round(bf$bf, 4), ".")
  })
  
  output$table_bf_linear <- renderDataTable({
    BF <- extractBF(linear_bf())
    
    bf_frame = get_table_content(BF$bf)
    
    colnames(bf_frame) <- c("Bayes factor", "interpretation")
    bf_frame
  }, rownames = FALSE, options = list(paging = FALSE, searching = FALSE))
  
})
