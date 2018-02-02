#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#
# App inspired by http://core-genomics.blogspot.com/2016/05/how-many-reads-to-sequence-genome.html

library(shiny)
library(dplyr)

# Define UI for application that draws a histogram
ui <- fluidPage(

   # Application title
   titlePanel("Coverage / Read Count Calculator"),

   h4("Calculate how much sequencing you need to hit a target depth of coverage (or vice versa)."),

   HTML("<p><strong>Instructions:</strong> set the read length/configuration and genome size, then select what you want to calculate.</p>"),

   HTML("<p>Written by <a href='http://stephenturner.us', target='blank'>Stephen Turner</a>, based on the <a href='http://www.ncbi.nlm.nih.gov/pubmed/3294162' target='_blank'>Lander-Waterman formula</a>, inspired by <a href='http://core-genomics.blogspot.com/2016/05/how-many-reads-to-sequence-genome.html' target='_blank'>a similar calculator</a> written by James Hadfield.<br/><em>C</em> = Coverage (X), <em>L</em> = Read length (bp), <em>G</em> = Haploid genome size (bp), <em>N</em> = Number of reads,  <em>ON</em> = On-target percentage and  <em>Dup</em> = Duplication percentage 
        <br/>Coverage is calculated as <em>C=LN/G*ON*(1-Dup)</em> and reads as <em>N=CG/L/ON/(1-Dup)</em>. <br/>Forked from <a href='https://github.com/stephenturner/covcalc' target='_blank'>GitHub</a><br/>
        Source code <a href='https://github.com/smw1414/covcalc' target='_blank'>GitHub</a></p>"),

   # HTML("<p>Written by <a href='http://stephenturner.us', target='blank'>Stephen Turner</a>, based on the <a href='http://www.ncbi.nlm.nih.gov/pubmed/3294162' target='_blank'>Lander-Waterman formula</a>, inspired by <a href='http://core-genomics.blogspot.com/2016/05/how-many-reads-to-sequence-genome.html' target='_blank'>a similar calculator</a> written by James Hadfield.</p><p>Coverage is calculated as <em>C=LN/G</em> and reads as <em>N=CG/L</em> where
   #        <ul>
   #          <li><em>C</em> = Coverage (X)</li>
   #          <li><em>L</em> = Read length (bp)</li>
   #          <li><em>G</em> = Haploid genome size (bp)</li>
   #          <li><em>N</em> = Number of reads</li>
   #        </ul>
   #      </p>"),

   # Sidebar with a slider input for number of bins
   sidebarLayout(
      sidebarPanel(
        fluidRow(column(
          6,
          numericInput(
            "readlength",
            label = h4("Read length (bp)"),
            value = 75
          ),
          numericInput(
            "ontarget",
            label = h4("Expected on-target (%)"),
            value = 70
          )
        ),
        column(
          6,
          numericInput(
            "dup",
            label = h4("Expected duplication (%)"),
            value = 15
          ),
          numericInput(
            "insert",
            label = h4("Average insert size (bp)"),
            value = 160
          )
        )), 
       
        radioButtons("sepe", label = NA,
                     choices = list("Paired-end" = 2, "Single-end" = 1),
                     selected = 2),



        HTML("<h4>Genome size</h4><p>Pick a genome below or manually enter the haploid genome size in bp. You can use scientific notation (e.g., enter <strong>3.2e9</strong> for 3,200,000,000bp, or 3.2 Gb).</p>"),

        radioButtons("genomesizebutton", label=NULL,
                     choices=list(
                       "Human (3.1 Gb)" = 3096649726,
                       "Agilent V6 exome (60 Mb)" = 60e6,
                       "TruSeq exome (45 Mb)" = 45e6,
                       "S. cerevisiae (12.2 Mb)" = 12157105,
                       "E. coli K-12 (4.6 Mb)" = 	4641652,
                       "Other (Enter manually)"=""),
                     selected=3096649726),

        uiOutput("genomesizeOtherUI")

      ),

      # Show a plot of the generated distribution
      mainPanel(

        wellPanel(
             radioButtons("whatcalc",
                          label=h3("What do you want to know?"),
                          choices=list("# Reads (how many reads do I need to hit a target depth of coverage?)"="reads",
                                       "Coverage (what's my coverage depth obtained from a set number of reads)"="coverage"),
                          selected="reads"),
             # hr(),
             # Dynamic input goes here
             uiOutput("ui")
        ),

        #h1(textOutput("text")),
        h1(htmlOutput("text")),
        
        # plotOutput("myplot"),

        tags$br()

      )
   )

)


# Define server logic required to draw a histogram
server <- function(input, output) {

  observe({
    if (input$genomesizebutton!="") {
      output$genomesizeOtherUI <- renderText(
        HTML("Selected genome size: <strong>", prettyNum(input$genomesizebutton, big.mark=","), "</strong>")
      )
    } else {
      output$genomesizeOtherUI <- renderUI({
        numericInput("genomesizeOther", label=NA, value=NULL)
      })
    }

  })

  genomesize <- reactive({
    if (input$genomesizebutton!="") {
      return(as.numeric(input$genomesizebutton))
    } else {
      return(input$genomesizeOther)
    }
  })

  # observe(print(genomesize()))

  output$ui <- renderUI({
    if (is.null(input$whatcalc)) {
      return()
    } else {
      switch(input$whatcalc,
             "reads"=sliderInput("dynamic", h4("Desired coverage"), min=5, max=500, value=30, step=5, post="x"),
             "coverage"=sliderInput("dynamic",
                                    label=h4("Number reads sequenced (millions)"),
                                    min=1, max=1001, value=384,
                                    post=paste0("M "))
      )
    }
  })

  # Converts genome size into SI-prefix character
  bptosi <- function(bp) {
    if (bp>1e9)      return(paste0(round(bp/1e9, 2), "Gb"))
    else if (bp>1e6) return(paste0(round(bp/1e6, 2), "Mb"))
    else if (bp>1e3) return(paste0(round(bp/1e3, 2), "Kb"))
    else             return(paste0(bp, "bp"))
  }

  outputtext <- reactive({
    req(input$whatcalc, genomesize)
    if (input$whatcalc=="reads") {
      paste0("Calculating reads", input$dynamic)
    } else if (input$whatcalc=="coverage") {
      paste0("Calculating coverage", input$dynamic)
    }
  })
  
  output$text <- renderText({
    #req(input$whatcalc)
    if (input$genomesizebutton=="") req(input$genomesizeOther)
    if (input$whatcalc=="reads") {
      # N=CG/L
      cycle <- as.numeric(input$readlength)*as.numeric(input$sepe)
      insert <- ifelse( cycle > input$insert, 
                       input$insert/cycle ,
                       1)
      reads <- (as.numeric(input$dynamic)*as.numeric(genomesize()))/cycle/as.numeric(input$ontarget*0.01)/as.numeric(1-(input$dup*0.01))/insert
      gb <- cycle*reads/1e+9
      str1 <- paste0(input$dynamic, "X coverage of ", bptosi(genomesize()), " target region")
      str2 <- paste0("Reads: ",signif(reads/1e6, 3)," million reads")
      str3 <- paste0("Read Length: ",input$readlength , " x ", input$sepe)
      str4 <- paste0("Bases output: ",round(gb,digits = 1),"Gb")
      HTML(paste(str1, str2, str3, str4, sep = '<br/>'))
      #reads <- (as.numeric(input$dynamic)*as.numeric(genomesize()))/(as.numeric(input$readlength)*as.numeric(input$sepe))       
      #paste0("*", signif(reads/1e6, 3), " million reads* required for ", input$dynamic, "X coverage of a ", bptosi(genomesize()), " genome using ", input$sepe, "x", input$readlength, " sequencing reads.")
    } else if (input$whatcalc=="coverage") {
      # C=LN/G
      cycle <- as.numeric(input$readlength)*as.numeric(input$sepe)
      insert <- ifelse( cycle > input$insert, 
                        input$insert/cycle ,
                        1)
      coverage <- (as.numeric(input$readlength)*as.numeric(input$sepe))*as.numeric(input$dynamic)*1e6/as.numeric(genomesize())*as.numeric(input$ontarget*0.01)*as.numeric(1-(input$dup*0.01))*insert
      gb <-  ((as.numeric(input$readlength)*as.numeric(input$sepe))*as.numeric(input$dynamic)*1e6)/1e9
      str1 <- paste0(round(coverage), "X coverage of ", bptosi(genomesize()), " target region")
      str2 <- paste0("Reads: ",input$dynamic," million reads")
      str3 <- paste0("Read Length: ",input$readlength , " x ", input$sepe)
      str4 <- paste0("Bases output: ",round(gb, digits = 1),"Gb")
      HTML(paste(str1, str2, str3, str4, sep = '<br/>'))
      #paste0("*", round(coverage), "X coverage* for a ", bptosi(genomesize()), " genome obtained with ", input$dynamic, "M ", input$sepe, "x", input$readlength, " sequencing reads.")
    }
  })

  # output$myplot <- renderPlot({
  #   req(input$dynamic)
  #   plot(input$dynamic)
  # })

   # output$distPlot <- renderPlot({
   #    # generate bins based on input$bins from ui.R
   #    x    <- faithful[, 2]
   #    bins <- seq(min(x), max(x), length.out = input$bins + 1)
   #
   #    # draw the histogram with the specified number of bins
   #    hist(x, breaks = bins, col = 'darkgray', border = 'white')
   # })

}

# Run the application
shinyApp(ui = ui, server = server)

