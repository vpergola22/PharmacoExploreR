# app.R
# Shiny app for PharmacoExploreR package
# Allows interactive exploration of pharmacogenomic data

library(shiny)
library(PharmacoExploreR)
library(PharmacoGx)
library(ggplot2)
library(DT)



# INCREASE FILE UPLOAD SIZE LIMIT

# Default is 5MB, increase to 500MB for PharmacoSet files
options(shiny.maxRequestSize = 500*1024^2)  # 500 MB


# UI

ui <- fluidPage(
  
  # Application title
  titlePanel("PharmacoExploreR: Interactive Pharmacogenomic Data Explorer"),
  
  # Add custom CSS
  tags$head(
    tags$style(HTML("
      .title-panel { background-color: #2c3e50; color: white; padding: 20px; }
      .well { background-color: #ecf0f1; }
      .btn-primary { background-color: #3498db; border-color: #2980b9; }
      .btn-primary:hover { background-color: #2980b9; }
    "))
  ),
  
  sidebarLayout(
    
    # create sidebar
    sidebarPanel(
      width = 3,
      
      h3("Data Input"),
      
      # Tab selection
      radioButtons(
        inputId = "dataSource",
        label = "Select Data Source:",
        choices = c(
          "Use Demo Data" = "demo",
          "Upload PharmacoSet (.rds)" = "upload"
        ),
        selected = "demo"
      ),
      
      # File upload (conditional)
      conditionalPanel(
        condition = "input.dataSource == 'upload'",
        fileInput(
          inputId = "psetFile",
          label = "Upload PharmacoSet (.rds file):",
          accept = c(".rds", ".RDS")
        ),
        helpText("Upload a PharmacoSet object saved as .rds file")
      ),
      
      # Demo data info (conditional)
      conditionalPanel(
        condition = "input.dataSource == 'demo'",
        helpText("Using included demo dataset with 15 drugs and 60 cell lines")
      ),
      
      hr(),
      
      # Drug selection
      uiOutput("drugSelector"),
      
      hr(),
      
      h3("Analysis Options"),
      
      # Correlation settings
      selectInput(
        inputId = "sensitivityMeasure",
        label = "Sensitivity Measure:",
        choices = c("AAC (Area Above Curve)" = "aac_recomputed",
                    "AUC (Area Under Curve)" = "auc_recomputed",
                    "IC50" = "ic50_recomputed"),
        selected = "aac_recomputed"
      ),
      
      selectInput(
        inputId = "correlationMethod",
        label = "Correlation Method:",
        choices = c("Pearson" = "pearson",
                    "Spearman" = "spearman",
                    "Kendall" = "kendall"),
        selected = "pearson"
      ),
      
      hr(),
      
      # Action button
      actionButton(
        inputId = "runAnalysis",
        label = "Run Analysis",
        class = "btn-primary btn-lg btn-block",
        icon = icon("play")
      ),
      
      hr(),
      
      # Download options
      downloadButton("downloadCorrelations", "Download Correlation Results"),
      br(), br(),
      downloadButton("downloadPlot", "Download Current Plot")
    ),
    
    # create main panel
    mainPanel(
      width = 9,
      
      tabsetPanel(
        id = "mainTabs",
        
        # TAB 1: Welcome
        tabPanel(
          "Welcome",
          icon = icon("home"),
          
          h2("Welcome to PharmacoExploreR!"),
          
          p("This interactive tool allows you to explore pharmacogenomic data 
            and identify relationships between gene expression and drug response."),
          
          h3("Getting Started"),
          
          tags$ol(
            tags$li(strong("Choose your data source:"), 
                    "Use the demo dataset or upload your own PharmacoSet (.rds file)"),
            tags$li(strong("Select a drug:"), 
                    "Choose from the dropdown menu in the sidebar"),
            tags$li(strong("Configure analysis:"), 
                    "Select sensitivity measure and correlation method"),
            tags$li(strong("Run analysis:"), 
                    "Click the 'Run Analysis' button"),
            tags$li(strong("Explore results:"), 
                    "Navigate through the tabs to view different visualizations")
          ),
          
          h3("Data Requirements"),
          
          p(strong("If uploading your own data:")),
          tags$ul(
            tags$li(strong("File type:"), ".rds file containing a PharmacoSet object"),
            tags$li(strong("Data structure:"), "Must be a valid PharmacoSet from PharmacoGx package"),
            tags$li(strong("Required data:"), 
                    "Must include molecular profiles (gene expression) and 
                    treatment response (drug sensitivity) data"),
            tags$li(strong("Example data:"), 
                    "Download demo data from ",
                    tags$a(href = "https://github.com/vpergola22/PharmacoExploreR/tree/main/data",
                           "GitHub repository"))
          ),
          
          h3("About the Demo Data"),
          
          p("The demo dataset is a synthetic pharmacogenomic dataset containing:"),
          tags$ul(
            tags$li("10 drugs with dose-response data"),
            tags$li("150 cell lines from 5 tissue types"),
            tags$li("Multiple sensitivity metrics (AAC, AUC, IC50)")
          ),
          
          h3("Analysis Features"),
          
          tags$ul(
            tags$li(strong("Correlation Analysis:"), 
                    "Identify genes associated with drug response"),
            tags$li(strong("Volcano Plots:"), 
                    "Visualize genome-wide gene-drug relationships"),
            tags$li(strong("Scatterplots:"), 
                    "Examine specific gene-drug correlations"),
            tags$li(strong("Response Groups:"), 
                    "Classify samples as sensitive or resistant"),
            tags$li(strong("Differential Expression:"), 
                    "Find genes differentially expressed between groups"),
            tags$li(strong("Dose-Response Curves:"), 
                    "Visualize drug effects across concentrations")
          ),
          
          hr(),
          
          p("Developed by Victoria Pergola for BCB410H: Applied Bioinformatics, 2025"),
          p("For more information, visit the ",
            tags$a(href = "https://github.com/vpergola22/PharmacoExploreR", 
                   "GitHub repository"))
        ),
        
        # TAB 2: Correlation Results
        tabPanel(
          "Correlation Results",
          icon = icon("table"),
          
          h3("Gene-Drug Correlation Analysis"),
          
          p("This table shows the correlation between each gene's expression 
            and drug sensitivity. Positive correlations indicate genes whose 
            higher expression is associated with drug sensitivity, while negative 
            correlations indicate resistance markers."),
          
          fluidRow(
            column(
              width = 12,
              DT::dataTableOutput("correlationTable")
            )
          ),
          
          hr(),
          
          fluidRow(
            column(
              width = 6,
              h4("Summary Statistics"),
              verbatimTextOutput("correlationSummary")
            ),
            column(
              width = 6,
              h4("Significant Genes"),
              verbatimTextOutput("significantGenes")
            )
          )
        ),
        
        # TAB 3: Volcano Plot
        tabPanel(
          "Volcano Plot",
          icon = icon("chart-area"),
          
          h3("Genome-Wide Correlation Overview"),
          
          p("The volcano plot displays correlation strength (x-axis) versus 
            statistical significance (y-axis) for all genes. Genes in the top 
            corners are both highly correlated and statistically significant."),
          
          fluidRow(
            column(
              width = 3,
              wellPanel(
                h4("Plot Settings"),
                numericInput(
                  "volcanoCorThreshold",
                  "Correlation Threshold:",
                  value = 0.3,
                  min = 0,
                  max = 1,
                  step = 0.05
                ),
                numericInput(
                  "volcanoPvalThreshold",
                  "P-value Threshold:",
                  value = 0.05,
                  min = 0.001,
                  max = 0.1,
                  step = 0.01
                )
              )
            ),
            column(
              width = 9,
              plotOutput("volcanoPlot", height = "600px")
            )
          )
        ),
        
        # TAB 4: Gene Scatterplot
        tabPanel(
          "Gene Scatterplot",
          icon = icon("chart-line"),
          
          h3("Gene Expression vs Drug Response"),
          
          p("Select a specific gene to visualize its expression relationship 
            with drug sensitivity across all cell lines."),
          
          fluidRow(
            column(
              width = 3,
              wellPanel(
                h4("Gene Selection"),
                uiOutput("geneSelector"),
                hr(),
                actionButton(
                  "selectTopGene",
                  "Select Top Correlated Gene",
                  class = "btn-primary btn-block"
                )
              )
            ),
            column(
              width = 9,
              plotOutput("scatterPlot", height = "600px")
            )
          )
        ),
        
        # TAB 5: Response Groups
        tabPanel(
          "Response Groups",
          icon = icon("object-group"),
          
          h3("Sample Classification"),
          
          p("Classify cell lines into sensitive and resistant groups based on 
            drug response."),
          
          fluidRow(
            column(
              width = 3,
              wellPanel(
                h4("Classification Method"),
                selectInput(
                  "groupMethod",
                  "Method:",
                  choices = c(
                    "Median Split" = "median",
                    "Quantile (Top/Bottom 25%)" = "quantile",
                    "Manual Threshold" = "manual"
                  ),
                  selected = "median"
                ),
                conditionalPanel(
                  condition = "input.groupMethod == 'manual'",
                  numericInput(
                    "manualThreshold",
                    "Threshold:",
                    value = 0.5,
                    min = 0,
                    max = 1,
                    step = 0.05
                  )
                ),
                hr(),
                actionButton(
                  "classifySamples",
                  "Classify Samples",
                  class = "btn-primary btn-block"
                )
              )
            ),
            column(
              width = 9,
              h4("Classification Results"),
              plotOutput("groupDistribution", height = "400px"),
              hr(),
              h4("Sample Group Assignments"),
              DT::dataTableOutput("groupTable")
            )
          )
        ),
        
        # TAB 6: Differential Expression
        tabPanel(
          "Differential Expression",
          icon = icon("exchange-alt"),
          
          h3("Differential Expression Analysis"),
          
          p("Identify genes differentially expressed between sensitive and 
            resistant cell lines."),
          
          fluidRow(
            column(
              width = 3,
              wellPanel(
                h4("Analysis Settings"),
                selectInput(
                  "diffExprMethod",
                  "Statistical Method:",
                  choices = c(
                    "T-test" = "t.test",
                    "Limma" = "limma"
                  ),
                  selected = "t.test"
                ),
                hr(),
                actionButton(
                  "runDiffExpr",
                  "Run Differential Expression",
                  class = "btn-primary btn-block"
                ),
                hr(),
                p(strong("Note:"), "You must classify samples first 
                  (see Response Groups tab)")
              )
            ),
            column(
              width = 9,
              h4("Top Differentially Expressed Genes"),
              DT::dataTableOutput("diffExprTable"),
              hr(),
              h4("Gene Expression by Group"),
              uiOutput("diffExprGeneSelector"),
              plotOutput("geneBoxplot", height = "500px")
            )
          )
        ),
        
        # TAB 7: Dose-Response
        tabPanel(
          "Dose-Response Curves",
          icon = icon("chart-line"),
          
          h3("Drug Dose-Response Curves"),
          
          p("Visualize how cell viability changes with drug concentration for 
            selected cell lines."),
          
          fluidRow(
            column(
              width = 3,
              wellPanel(
                h4("Cell Line Selection"),
                uiOutput("cellLineSelector"),
                hr(),
                actionButton(
                  "plotDoseResponse",
                  "Plot Dose-Response",
                  class = "btn-primary btn-block"
                )
              )
            ),
            column(
              width = 9,
              plotOutput("doseResponsePlot", height = "600px")
            )
          )
        ),
        
        # TAB 8: Help
        tabPanel(
          "Help",
          icon = icon("question-circle"),
          
          h2("Help & Documentation"),
          
          h3("Interpreting Results"),
          
          h4("Correlation Analysis"),
          tags$ul(
            tags$li(strong("Correlation coefficient (r):"), 
                    "Ranges from -1 to 1. Values close to 0 indicate weak 
                    correlation, while values close to ±1 indicate strong correlation."),
            tags$li(strong("P-value:"), 
                    "Indicates statistical significance. Values < 0.05 are 
                    typically considered significant."),
            tags$li(strong("Adjusted P-value:"), 
                    "P-values corrected for multiple testing using FDR 
                    (False Discovery Rate).")
          ),
          
          h4("Sensitivity Measures"),
          tags$ul(
            tags$li(strong("AAC (Area Above Curve):"), 
                    "Higher values indicate greater drug sensitivity."),
            tags$li(strong("AUC (Area Under Curve):"), 
                    "Lower values indicate greater drug sensitivity."),
            tags$li(strong("IC50:"), 
                    "Drug concentration required to inhibit growth by 50%. 
                    Lower values indicate greater sensitivity.")
          ),
          
          h3("Common Issues"),
          
          h4("No correlations showing?"),
          tags$ul(
            tags$li("Ensure you've clicked 'Run Analysis'"),
            tags$li("Check that the selected drug has data in your PharmacoSet"),
            tags$li("Try a different drug or sensitivity measure")
          ),
          
          h4("Upload fails?"),
          tags$ul(
            tags$li("Ensure file is a valid PharmacoSet saved as .rds"),
            tags$li("File must contain both molecular profiles and 
                    sensitivity data"),
            tags$li("Maximum file size: 100 MB")
          ),
          
          h4("Plots not displaying?"),
          tags$ul(
            tags$li("Run the analysis first using 'Run Analysis' button"),
            tags$li("For some plots, additional steps may be required 
                    (e.g., classify samples first)")
          ),
          
          h3("Getting Example Data"),
          
          p("Demo data is included with the app. To download or explore:"),
          tags$ul(
            tags$li("Visit: ", 
                    tags$a(href = "https://github.com/vpergola22/PharmacoExploreR/tree/main/data",
                           "GitHub Repository")),
            tags$li("File: nci60_subset.rds"),
            tags$li("Or use PharmacoGx to download public datasets")
          ),
          
          h3("References"),
          
          p("Smirnov, P., et al. (2016). PharmacoGx: an R package for 
            analysis of large pharmacogenomic datasets. Bioinformatics, 32(8), 1244-1246."),
          
          p("For package documentation, visit: ",
            tags$a(href = "https://github.com/vpergola22/PharmacoExploreR",
                   "https://github.com/vpergola22/PharmacoExploreR"))
        )
      )
    )
  )
)


# SERVER

server <- function(input, output, session) {
  
  # REACTIVE VALUES
  
  rv <- reactiveValues(
    pset = NULL,
    correlations = NULL,
    responseGroups = NULL,
    diffExprResults = NULL
  )
  
  # LOAD DATA
  
  # Load PharmacoSet (demo or uploaded)
  observe({
    if (input$dataSource == "demo") {
      # Load demo data
      tryCatch({
        rv$pset <- readRDS(system.file("extdata", "nci60_subset.rds",
                                       package = "PharmacoExploreR"))
        showNotification("Demo data loaded successfully!", 
                         type = "message", duration = 3)
      }, error = function(e) {
        showNotification(paste("Error loading demo data:", e$message),
                         type = "error", duration = 10)
      })
    }
  })
  
  # Handle file upload
  observeEvent(input$psetFile, {
    req(input$psetFile)
    
    tryCatch({
      rv$pset <- readRDS(input$psetFile$datapath)
      
      # Validate it's a PharmacoSet
      if (!inherits(rv$pset, "PharmacoSet")) {
        showNotification("Uploaded file is not a valid PharmacoSet!",
                         type = "error", duration = 10)
        rv$pset <- NULL
      } else {
        showNotification("PharmacoSet uploaded successfully!",
                         type = "message", duration = 3)
      }
    }, error = function(e) {
      showNotification(paste("Error uploading file:", e$message),
                       type = "error", duration = 10)
    })
  })
  
  # DYNAMIC UI ELEMENTS
  
  # Drug selector
  output$drugSelector <- renderUI({
    req(rv$pset)
    
    drugs <- drugNames(rv$pset)
    
    selectInput(
      inputId = "selectedDrug",
      label = "Select Drug:",
      choices = drugs,
      selected = drugs[1]
    )
  })
  
  # Gene selector
  output$geneSelector <- renderUI({
    req(rv$correlations)
    
    genes <- rv$correlations$gene
    
    selectInput(
      inputId = "selectedGene",
      label = "Select Gene:",
      choices = genes,
      selected = genes[1]
    )
  })
  
  # Diff expr gene selector
  output$diffExprGeneSelector <- renderUI({
    req(rv$diffExprResults)
    
    genes <- rv$diffExprResults$gene
    
    selectInput(
      inputId = "diffExprSelectedGene",
      label = "Select Gene to Visualize:",
      choices = genes,
      selected = genes[1]
    )
  })
  
  # Cell line selector
  output$cellLineSelector <- renderUI({
    req(rv$pset)
    
    cellLines <- cellNames(rv$pset)
    
    checkboxGroupInput(
      inputId = "selectedCellLines",
      label = "Select Cell Lines (max 5):",
      choices = cellLines,
      selected = cellLines[1:min(3, length(cellLines))]
    )
  })
  
  # RUN ANALYSIS
  
  observeEvent(input$runAnalysis, {
    req(rv$pset, input$selectedDrug)
    
    withProgress(message = 'Running analysis...', value = 0, {
      
      tryCatch({
        incProgress(0.3, detail = "Computing correlations...")
        
        rv$correlations <- correlateExpressionAUC(
          pset = rv$pset,
          drug = input$selectedDrug,
          mDataType = "rna",
          sensitivity.measure = input$sensitivityMeasure,
          method = input$correlationMethod
        )
        
        incProgress(0.7, detail = "Analysis complete!")
        
        showNotification("Analysis completed successfully!",
                         type = "message", duration = 3)
        
        # Switch to correlation results tab
        updateTabsetPanel(session, "mainTabs", selected = "Correlation Results")
        
      }, error = function(e) {
        showNotification(paste("Error in analysis:", e$message),
                         type = "error", duration = 10)
      })
    })
  })
  
  # CORRELATION TABLE
  
  output$correlationTable <- DT::renderDataTable({
    req(rv$correlations)
    
    DT::datatable(
      rv$correlations,
      options = list(
        pageLength = 15,
        scrollX = TRUE,
        order = list(list(2, 'desc'))  # Sort by correlation
      ),
      rownames = FALSE
    ) %>%
      DT::formatRound(columns = c('cor', 'pval', 'adj_pval'), digits = 4) %>%
      DT::formatStyle(
        'cor',
        backgroundColor = styleInterval(
          cuts = c(-0.3, 0.3),
          values = c('#ffcccc', 'white', '#ccffcc')
        )
      )
  })
  
  # Correlation summary
  output$correlationSummary <- renderPrint({
    req(rv$correlations)
    
    cat("Total genes analyzed:", nrow(rv$correlations), "\n")
    cat("Correlation range:", 
        round(range(rv$correlations$cor, na.rm = TRUE), 3), "\n")
    cat("Mean |correlation|:",
        round(mean(abs(rv$correlations$cor), na.rm = TRUE), 3), "\n")
  })
  
  # Significant genes
  output$significantGenes <- renderPrint({
    req(rv$correlations)
    
    sig <- rv$correlations[rv$correlations$adj_pval < 0.05, ]
    
    cat("Significant genes (adj p < 0.05):", nrow(sig), "\n")
    cat("Positively correlated:", 
        sum(sig$cor > 0, na.rm = TRUE), "\n")
    cat("Negatively correlated:",
        sum(sig$cor < 0, na.rm = TRUE), "\n")
  })
  
  # VOLCANO PLOT
  
  output$volcanoPlot <- renderPlot({
    req(rv$correlations)
    
    volcanoAUC(
      corResults = rv$correlations,
      corThreshold = input$volcanoCorThreshold,
      pvalThreshold = input$volcanoPvalThreshold
    )
  })
  
  # SCATTERPLOT
  
  # Select top gene button
  observeEvent(input$selectTopGene, {
    req(rv$correlations)
    
    top_gene <- rv$correlations$gene[which.max(abs(rv$correlations$cor))]
    updateSelectInput(session, "selectedGene", selected = top_gene)
  })
  
  output$scatterPlot <- renderPlot({
    req(rv$pset, rv$correlations, input$selectedGene, input$selectedDrug)
    
    plotExprAUC(
      pset = rv$pset,
      corResults = rv$correlations,
      gene = input$selectedGene,
      drug = input$selectedDrug,
      mDataType = "rna",
      sensitivity.measure = input$sensitivityMeasure,
      method = input$correlationMethod
    )
  })
  
  # RESPONSE GROUPS
  
  observeEvent(input$classifySamples, {
    req(rv$pset, input$selectedDrug)
    
    withProgress(message = 'Classifying samples...', value = 0.5, {
      
      tryCatch({
        threshold_val <- if (input$groupMethod == "manual") {
          input$manualThreshold
        } else {
          NULL
        }
        
        rv$responseGroups <- defineResponseGroups(
          pset = rv$pset,
          drug = input$selectedDrug,
          sensitivity.measure = input$sensitivityMeasure,
          method = input$groupMethod,
          threshold = threshold_val
        )
        
        showNotification("Samples classified successfully!",
                         type = "message", duration = 3)
        
      }, error = function(e) {
        showNotification(paste("Error classifying samples:", e$message),
                         type = "error", duration = 10)
      })
    })
  })
  
  output$groupDistribution <- renderPlot({
    req(rv$responseGroups)
    
    df <- data.frame(Group = rv$responseGroups)
    
    ggplot(df, aes(x = Group, fill = Group)) +
      geom_bar() +
      scale_fill_manual(values = c("sensitive" = "#2ecc71", 
                                   "resistant" = "#e74c3c")) +
      labs(title = "Sample Distribution by Response Group",
           x = "Response Group",
           y = "Number of Samples") +
      theme_minimal(base_size = 14) +
      theme(legend.position = "none")
  })
  
  output$groupTable <- DT::renderDataTable({
    req(rv$responseGroups)
    
    df <- data.frame(
      Sample = names(rv$responseGroups),
      Group = as.character(rv$responseGroups)
    )
    
    DT::datatable(
      df,
      options = list(pageLength = 10, scrollX = TRUE),
      rownames = FALSE
    ) %>%
      DT::formatStyle(
        'Group',
        backgroundColor = styleEqual(
          c('sensitive', 'resistant'),
          c('#d5f4e6', '#fadbd8')
        )
      )
  })
  
  # DIFFERENTIAL EXPRESSION
  
  observeEvent(input$runDiffExpr, {
    req(rv$pset, rv$responseGroups)
    
    withProgress(message = 'Running differential expression...', value = 0.5, {
      
      tryCatch({
        rv$diffExprResults <- runDiffExpr(
          pset = rv$pset,
          groupLabels = rv$responseGroups,
          mDataType = "rna",
          method = input$diffExprMethod
        )
        
        showNotification("Differential expression completed!",
                         type = "message", duration = 3)
        
      }, error = function(e) {
        showNotification(paste("Error in differential expression:", e$message),
                         type = "error", duration = 10)
      })
    })
  })
  
  output$diffExprTable <- DT::renderDataTable({
    req(rv$diffExprResults)
    
    # Sort by adj_pval
    df_sorted <- rv$diffExprResults[order(rv$diffExprResults$adj_pval), ]
    
    DT::datatable(
      head(df_sorted, 50),  # Top 50 genes
      options = list(pageLength = 15, scrollX = TRUE),
      rownames = FALSE
    ) %>%
      DT::formatRound(columns = c('logFC', 'pval', 'adj_pval'), digits = 4) %>%
      DT::formatStyle(
        'logFC',
        backgroundColor = styleInterval(
          cuts = c(-1, 1),
          values = c('#ffcccc', 'white', '#ccffcc')
        )
      )
  })
  
  output$geneBoxplot <- renderPlot({
    req(rv$pset, rv$responseGroups, input$diffExprSelectedGene)
    
    plotGeneBoxplot(
      pset = rv$pset,
      gene = input$diffExprSelectedGene,
      groupLabels = rv$responseGroups,
      mDataType = "rna",
      plotType = "boxplot"
    )
  })
  
  # DOSE-RESPONSE
  
  output$doseResponsePlot <- renderPlot({
    req(rv$pset, input$selectedCellLines)
    
    # Limit to 5 cell lines
    selected <- input$selectedCellLines[1:min(5, length(input$selectedCellLines))]
    
    plotDoseResponse(
      pset = rv$pset,
      cell.lines = selected,
      sensitivity.measure = "Viability"
    )
  })
  
  # DOWNLOAD HANDLERS
  
  output$downloadCorrelations <- downloadHandler(
    filename = function() {
      paste0("correlations_", input$selectedDrug, "_", 
             Sys.Date(), ".csv")
    },
    content = function(file) {
      req(rv$correlations)
      write.csv(rv$correlations, file, row.names = FALSE)
    }
  )
  
  output$downloadPlot <- downloadHandler(
    filename = function() {
      paste0("plot_", Sys.Date(), ".png")
    },
    content = function(file) {
      # Save the current plot
      # This is a placeholder - implement based on current tab
      showNotification("Plot download feature coming soon!",
                       type = "message")
    }
  )
}


# RUN APP

shinyApp(ui = ui, server = server)