library(shiny)
library(plotly)
library(bslib)
library(DBI)
library(RSQLite)
library(DT)
source("scripts/plot_section.R")
source("scripts/search_section.R")

ui <- page_navbar(
  id    = "main_navbar",
  title = "mBAT-combo Portal",
  theme = bs_theme(
    bootswatch = "cosmo",
    primary = "#007BFF"
  ),
  
  bg = "#E9EEF5",
  fillable = FALSE,
  
  header = tags$style(HTML("
  /* Whole navbar */
  .navbar {
    background-color: #E9EEF5 !important;
    min-height: 92px;
    padding: 0 28px !important;
    border-bottom: 8px solid #073B78;
    box-shadow: none;
  }

  /* App title */
  .navbar-brand {
    color: #073B78 !important;
    font-size: 2.4rem;
    font-weight: 800;
    letter-spacing: -0.04em;
    margin-right: 80px;
  }

  /* Navbar layout */
  .navbar-nav {
    align-items: stretch;
    gap: 42px;
  }

  /* Tabs: Home, Results, Search */
  .navbar-nav .nav-link {
    position: relative;
    color: #073B78 !important;
    font-size: 1.15rem;
    font-weight: 700;
    padding: 32px 4px 28px 4px !important;
    background: transparent !important;
    border-radius: 0 !important;
    box-shadow: none !important;
  }

  .navbar-nav .nav-link.active {
    background: transparent !important;
    color: #073B78 !important;
  }

  /* Hover */
  .navbar-nav .nav-link:hover {
    color: #007BFF !important;
    background: transparent !important;
  }
")),
  
  nav_panel("Home",
            div(
              style = paste0(
                "position: relative; min-height: calc(100vh - 56px);",
                "background: url('main.png') center/cover no-repeat #0d2d6b;",
                "display: flex; flex-direction: column; justify-content: center; padding: 0 8%;"
              ),
              div(
                style = "max-width: 680px;",
                tags$h1(
                  "mBAT-combo Portal",
                  style = "color: #ffffff; font-size: 2.4rem; font-weight: 700; margin-bottom: 1rem;"
                ),
                tags$p(
                  paste0(
                    "mBAT-combo Portal facilitates integrative gene-based association analysis ",
                    "of GWAS summary statistics. It applies the mBAT-combo, mBAT, and fastBAT ",
                    "methods to identify genes associated with complex traits, and provides ",
                    "interactive exploration of results across thousands of UK Biobank phenotypes."
                  ),
                  style = "color: rgba(255,255,255,0.85); font-size: 1.1rem; line-height: 1.7; margin-bottom: 2rem;"
                ),
                div(
                  style = "display: flex; gap: 12px; flex-wrap: wrap;",
                  actionButton("home_browse", "Browse Results",
                               class = "btn btn-primary btn-lg",
                               style = "background:#2c6fad; border-color:#2c6fad;"),
                  actionButton("home_search", "Search Genes / Traits",
                               class = "btn btn-outline-light btn-lg")
                )
              )
            )
  ),
  
  nav_panel("Results",
              uiOutput("trait_title"),
              card(card_header("Manhattan Plot — P_mBATcombo"), plotlyOutput("manhattan_mBATcombo")),
              card(card_header("Manhattan Plot — P_mBAT"),      plotlyOutput("manhattan_mBAT")),
              card(card_header("Manhattan Plot — P_fastBAT"),   plotlyOutput("manhattan_fastBAT")),
              card(
                card_header(
                  "Top Genes",
                  downloadButton("download_csv", "Download CSV", class = "btn-sm btn-outline-primary float-end")
                ),
                DTOutput("gene_table")
              )
  ),
  
  nav_panel("Search",
            div(class = "container mt-5",
                h1("Search"),
                div(class = "input-group input-group-lg",
                    tags$select(
                      id       = "search_type",
                      class    = "form-select",
                      style    = "max-width:130px;",
                      onchange = "document.getElementById('filter_btn').style.display = this.value === 'gene' ? '' : 'none';",
                      tags$option(value = "gene",  "Gene"),
                      tags$option(value = "trait", "Trait")
                    ),
                    tags$input(
                      id          = "gene_query",
                      type        = "text",
                      class       = "form-control",
                      placeholder = "Search gene symbol or name..."
                    ),
                    actionButton("filter_btn", label = tagList(icon("sliders"), " Filters"),
                                 class = "btn-outline-secondary")
                ),
                # to check for the enter from the search box. 
                tags$script(HTML("document.getElementById('gene_query').addEventListener('keydown', function(e) 
                {if (e.key === 'Enter') {e.preventDefault();
                Shiny.setInputValue('search_submit', {query: this.value,nonce: Math.random()}, {priority: 'event'});}});")),
                uiOutput("filter_panel"),
                uiOutput("gene_suggestions"), 
                uiOutput("trait_suggestions"),
                uiOutput("search_results")
            )
  )
)

server <- function (input, output, session){
  
  #------------------------------
  # the section for graph
  #------------------------------
  selected_trait <- reactiveVal("100001")
  
  output$trait_title <- renderUI({
    tid   <- selected_trait()
    trait_meta <- readRDS("data/trait_metadata.rds")
    tname <- trait_meta$trait_name[trait_meta$trait_id == tid]
    rm(trait_meta)
    if (length(tname) == 0 || is.na(tname)) tname <- tid
    div(
      h2(tname),
      tags$span(tid, style = "color:#6c757d; font-size:0.9em;")
    )
  })
  
  output$manhattan_mBATcombo <- renderPlotly({
    plot(file.path("data", paste0(selected_trait(), ".rds")), "P_mBATcombo")
  })
  
  output$manhattan_mBAT <- renderPlotly({
    plot(file.path("data", paste0(selected_trait(), ".rds")), "P_mBAT")
  })
  
  output$manhattan_fastBAT <- renderPlotly({
    plot(file.path("data", paste0(selected_trait(), ".rds")), "P_fastBAT")
  })
  
  output$gene_table <- renderDT({
    data <- readRDS(file.path("data", paste0(selected_trait(), ".rds")))
    data[order(data$P_mBATcombo), c("Symbol", "Chr", "Start", "End", "No.SNPs",
                                    "P_mBATcombo", "P_mBAT", "P_fastBAT",
                                    "TopSNP", "TopSNP_Pvalue", "No.Eigenvalues")]
    datatable(data, rownames = FALSE,
              options  = list(
                searching    = FALSE,
                lengthChange = FALSE,
                pageLength   = 25,
                scrollX      = TRUE,
                scrollCollapse = TRUE
              ))
  })
  
  output$download_csv <- downloadHandler(
    filename = function() paste0(selected_trait(), ".csv"),
    content  = function(download_file) {
      data <- readRDS(file.path("data", paste0(selected_trait(), ".rds")))
      write.csv(data, download_file, row.names = FALSE)
    })
  
  #---------------------------
  # section for search
  #---------------------------
  
  
  gene_query_d <- debounce(reactive(input$gene_query), 300) 
  filter_open   <- reactiveVal(FALSE)
  gene_search_data <- reactiveVal(NULL)
  trait_search_data <- reactiveVal(NULL) 
  
  
  observeEvent(input$filter_btn, {
    filter_open(!filter_open())
  })
  
  output$filter_panel <- renderUI({
    
    # To avoid junk UI. 
    if (!filter_open()) return(NULL)
    
    div(
      class = "card mt-1 p-3",
      style = "max-width:360px; margin-left:auto;", 
      filter_row("P_mBATcombo", "cb_p1", "thresh_p1"),
      filter_row("P_mBAT",      "cb_p2", "thresh_p2"),
      filter_row("P_fastBAT",   "cb_p3", "thresh_p3")
    )
  })
  
  
  # Two functions for rendering the UI for the drop down suggestions. 
  output$gene_suggestions <- renderUI({
    
    # To avoid junk UI. 
    if (filter_open()) return(NULL)
    
    if (is.null(input$search_type) || input$search_type == "trait") return(NULL)
    
    q <- gene_query_d()
    if (nchar(q) < 2) {
      return(NULL)
    }

    div(
      style = "border:1px solid #dee2e6; border-radius:4px; max-height:260px; overflow-y:auto; background:white;",
      search_results_rows(search_gene(q), "gene")
    )
    
  })
  
  output$trait_suggestions <- renderUI({
    
    if (is.null(input$search_type) || input$search_type != "trait") return(NULL)
    
    q <- gene_query_d()
    
    if (nchar(q) < 2) {
      return(NULL)
    }
    
    div(
      style = "border:1px solid #dee2e6; border-radius:4px; max-height:260px; overflow-y:auto; background:white;",
      # Actually the tag does not matter can be aaaa. 
      search_results_rows(search_trait(q), "trait")
    )
    
  })
  
  # ------------- What happens if the user press enter. Handles for both trait and gene ---------------
  observeEvent(input$search_submit, {
    
    q <- input$search_submit$query
    if (is.null(q)) q <- ""
    q <- trimws(q)
    
    if (nchar(q) < 2) return()
    
    updateTextInput(session, "gene_query", value = "")
    
    search_type <- input$search_type
    if (is.null(search_type)) search_type <- "gene"
    
    if (search_type == "trait") {
      traits <- search_trait(q)
      
      trait_search_data(traits)
      gene_search_data(NULL)
      
      return()
    }
    
    genes <- search_gene(q)
    if (nrow(genes) == 0) {
      gene_search_data(data.frame())
      trait_search_data(NULL)
      return()
    }
    
    p1 <- if (isTRUE(input$cb_p1)) input$thresh_p1 else 1
    p2 <- if (isTRUE(input$cb_p2)) input$thresh_p2 else 1
    p3 <- if (isTRUE(input$cb_p3)) input$thresh_p3 else 1
    
    results <- reverse_index(genes$id, p1, p2, p3)
    
    gene_search_data(results)
    trait_search_data(NULL)
  }, ignoreInit = TRUE)
  
  
  # -------- What if a drop down suggestion result is clicked. Two functions for trait and gene. -----------------------
  
  observeEvent(input$gene_clicked, {
    
    gene_id <- input$gene_clicked
    if (is.null(gene_id) || !nzchar(gene_id)) return()
    
    updateTextInput(session, "gene_query", value = "")
    
    p1 <- if (isTRUE(input$cb_p1)) input$thresh_p1 else 1
    p2 <- if (isTRUE(input$cb_p2)) input$thresh_p2 else 1
    p3 <- if (isTRUE(input$cb_p3)) input$thresh_p3 else 1
    
    results <- reverse_index(gene_id, p1, p2, p3)
    gene_search_data(results)
  }, ignoreInit = TRUE)
  
  observeEvent(input$trait_click, {
    trait_search_data(NULL)
    gene_search_data(NULL)
    
    selected_trait(input$trait_click)
    updateNavbarPage(session, "main_navbar", selected = "Results")
  })
  
  # ----------------------- code that renders the output results -------------
  
  output$search_results <- renderUI({
    if (input$search_type == "gene"){
      render_results(gene_search_data(), "gene")
    } else if (input$search_type == "trait"){
      render_results(trait_search_data(), "trait")
    }
  })
  
  # ----------------
  # home page section
  #-----------------
  
  observeEvent(input$home_browse, {
    updateNavbarPage(session, "main_navbar", selected = "Results")
  })
  observeEvent(input$home_search, {
    updateNavbarPage(session, "main_navbar", selected = "Search")
  })
  
  observeEvent(input$trait_click, {
    selected_trait(input$trait_click)
    updateNavbarPage(session, "main_navbar", selected = "Results")
  })
  
}

shinyApp(ui = ui, server = server)