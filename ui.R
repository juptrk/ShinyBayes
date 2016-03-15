library(shiny)

# ui.R
shinyUI(navbarPage("Bayes introduction",
                   tabPanel("Start",
                            includeMarkdown("start.Rmd")),
                   tabPanel("Frequentist vs. Bayesian",
                            includeMarkdown("information.Rmd")),
                   tabPanel("Coin Flip Data",
                            sidebarLayout(
                              sidebarPanel(
                                "Our Beta distribution is beta(1,1).",
                                sliderInput("k",
                                            label = "k",
                                            min = 1,
                                            max = 10,
                                            value = 7),
                                sliderInput("N",
                                            label = "N",
                                            min = 1,
                                            max = 30,
                                            value = 25),
                                sliderInput("theta",
                                            label = "theta",
                                            min = 0,
                                            max = 1,
                                            value = .5,
                                            step = .025),
                                h4("estimation"),
                                sliderInput("hdi",
                                            label = "HDI",
                                            min = 50,
                                            max = 99.5,
                                            step = .5,
                                            value = 95),
                                numericInput("ROPE",
                                             label = "ROPE epsilon",
                                             min = 0,
                                             max = .5,
                                             value = .1),
                                h4("p-value"),
                                selectInput("option",
                                            label = "options",
                                            choices = c("k fix", "N fix"))
                              ),
                              mainPanel(
                                tabsetPanel(type = "tabs", 
                                            tabPanel("p-value",
                                                     textOutput("text_pvalue"),
                                                     plotOutput("plot_pvalue"),
                                                     includeMarkdown("pvalue.Rmd")),
                                            tabPanel("comparison",
                                                     textOutput("text_comparison"),
                                                     dataTableOutput("table_bf"),
                                                     includeMarkdown("comparison.Rmd")),
                                            tabPanel("estimation",
                                                     plotOutput("plot_estimation"),
                                                     textOutput("text_rope"),
                                                     includeMarkdown("estimation.Rmd"))
                                )
                              )
                            )
                   ),
                   tabPanel("Bivariate Normal Data",
                            sidebarLayout(
                              sidebarPanel(
                                sliderInput("N_cor",
                                            label = "Specify the number of samples",
                                            min = 1,
                                            max = 150,
                                            value = 50),
                                sliderInput("corr",
                                            label = "Specify the correlation",
                                            min = 0,
                                            max = 1,
                                            value = .5,
                                            step = .01),
                                textInput("Xaxis",
                                          label = "Set the text for X and Y axis",
                                          value = "X"),
                                textInput("Yaxis",
                                          label = "",
                                          value = "Y")
                              ),
                              mainPanel(
                                tabsetPanel(type = "tabs",
                                            tabPanel("correlation",
                                                     plotOutput("cor_plot")
                                            ),
                                            tabPanel("Bayes factor",
                                                     p("As we also explain in the tab Frequentist vs. Bayesian,
                                                     we calculate the Bayes factor for deciding, whether the
                                                     data fits to a linear model or not."),
                                                     p("We assume for our linear model, that the Event Y is
                                                      linear related to the Event X."),
                                                     p("In the tabular you can see, what the interpration is
                                                       for the calculated Bayes factor."),
                                                     p("A big Bayes factor confirmes the linear model in this case."),
                                                     textOutput("linear_bayes_text"),
                                                     dataTableOutput("table_bf_linear"))
                                )
                              )
                              
                            )
                   )
))

