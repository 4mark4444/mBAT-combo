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
                conditionalPanel(
                  condition = "output.filter_visible == true",
                  div(
                    class = "card mt-1 p-3",
                    style = "max-width:360px; margin-left:auto;",
                    filter_row("P_mBATcombo", "cb_p1", "thresh_p1"),
                    filter_row("P_mBAT",      "cb_p2", "thresh_p2"),
                    filter_row("P_fastBAT",   "cb_p3", "thresh_p3")
                  )
                ),
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
  selected_trait <- reactiveVal("Q_100001")
  highlight_gene <- reactiveVal(NULL)   # ensembl id to star on the plots, or NULL
  
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
    plot(file.path("data", paste0(selected_trait(), ".rds")), "P_mBATcombo", highlight_gene())
  })

  output$manhattan_mBAT <- renderPlotly({
    plot(file.path("data", paste0(selected_trait(), ".rds")), "P_mBAT", highlight_gene())
  })

  output$manhattan_fastBAT <- renderPlotly({
    plot(file.path("data", paste0(selected_trait(), ".rds")), "P_fastBAT", highlight_gene())
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
  filter_open  <- reactiveVal(FALSE)

  # ---- gene-search state (two-level lazy: paged cards + per-card load-more) ----
  gene_hits    <- reactiveVal(NULL)            # data.frame(id, symbol, genename)
  gene_page    <- reactiveVal(1)               # current page of gene cards
  gene_filters <- reactiveVal(list(p1 = NA, p2 = NA, p3 = NA))
  gstate       <- reactiveVal(list())          # gene id -> list(rows, cursor, more)

  # ---- trait-search state (render pagination only) ----
  trait_hits   <- reactiveVal(NULL)            # data.frame(id, trait_name)
  trait_page   <- reactiveVal(1)

  # Current filter thresholds; NA means that filter is inactive (checkbox unticked).
  current_filters <- function() {
    active <- function(cb, num) {
      isTRUE(cb) && length(num) && !is.na(num)
    }
    list(
      p1 = if (active(input$cb_p1, input$thresh_p1)) input$thresh_p1 else NA,
      p2 = if (active(input$cb_p2, input$thresh_p2)) input$thresh_p2 else NA,
      p3 = if (active(input$cb_p3, input$thresh_p3)) input$thresh_p3 else NA
    )
  }

  # Lazily fetch the first page of traits for any of `ids` not already loaded.
  ensure_gene_loaded <- function(ids) {
    s <- gstate(); f <- gene_filters(); changed <- FALSE
    for (id in ids) {
      if (is.null(s[[id]])) {
        res  <- reverse_index_page(id, f$p1, f$p2, f$p3, after_rowid = 0, limit = PAGE_SIZE + 1)
        more <- nrow(res) > PAGE_SIZE
        if (more) res <- res[seq_len(PAGE_SIZE), , drop = FALSE]
        s[[id]] <- list(rows   = res,
                        cursor = if (nrow(res)) res$rid[nrow(res)] else 0,
                        more   = more)
        changed <- TRUE
      }
    }
    if (changed) gstate(s)
  }

  observeEvent(input$filter_btn, {
    filter_open(!filter_open())
  })
  
  # The filter panel itself lives statically in the UI (inside a conditionalPanel)
  # so its checkbox/number inputs are created ONCE and keep their state across
  # open/close. This output just drives that panel's visibility.
  output$filter_visible <- reactive({ filter_open() })
  outputOptions(output, "filter_visible", suspendWhenHidden = FALSE)
  
  
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
      trait_hits(search_trait(q))
      trait_page(1)
      gene_hits(NULL)
      return()
    }

    # gene mode: resolve matches, reset state, eagerly load only page 1's genes.
    genes <- search_gene(q)
    gene_filters(current_filters())
    gstate(list())                 # drop any per-gene state from a previous search
    gene_page(1)
    gene_hits(genes)
    trait_hits(NULL)
    if (nrow(genes) > 0) {
      ensure_gene_loaded(genes$id[seq_len(min(PAGE_SIZE, nrow(genes)))])
    }
  }, ignoreInit = TRUE)
  
  
  # -------- What if a drop down suggestion result is clicked. Two functions for trait and gene. -----------------------
  
  observeEvent(input$gene_clicked, {
    
    gene_id <- input$gene_clicked
    if (is.null(gene_id) || !nzchar(gene_id)) return()
    
    updateTextInput(session, "gene_query", value = "")
    
    gene_filters(current_filters())
    gstate(list())
    gene_page(1)
    gene_hits(gene_row(gene_id))
    trait_hits(NULL)
    ensure_gene_loaded(gene_id)
    
  }, ignoreInit = TRUE)
  
  # ---- load 10 more traits for one gene card (rowid cursor, appends) ----
  observeEvent(input$load_more, {
    id <- input$load_more$id
    if (is.null(id)) return()
    s  <- gstate(); st <- s[[id]]
    if (is.null(st) || !isTRUE(st$more)) return()
    f   <- gene_filters()
    res <- reverse_index_page(id, f$p1, f$p2, f$p3,
                              after_rowid = st$cursor, limit = PAGE_SIZE + 1)
    more <- nrow(res) > PAGE_SIZE
    if (more) res <- res[seq_len(PAGE_SIZE), , drop = FALSE]
    s[[id]] <- list(rows   = rbind(st$rows, res),
                    cursor = if (nrow(res)) res$rid[nrow(res)] else st$cursor,
                    more   = more)
    gstate(s)
  }, ignoreInit = TRUE)

  # ---- page through gene cards (lazily load the new page's genes) ----
  observeEvent(input$gene_page, {
    hits <- gene_hits()
    if (is.null(hits) || nrow(hits) == 0) return()
    npages <- max(1, ceiling(nrow(hits) / PAGE_SIZE))
    n <- max(1, min(as.integer(input$gene_page$n), npages))
    gene_page(n)
    idx <- ((n - 1) * PAGE_SIZE + 1):min(n * PAGE_SIZE, nrow(hits))
    ensure_gene_loaded(hits$id[idx])
  }, ignoreInit = TRUE)

  # ---- page through trait cards (pure render pagination) ----
  observeEvent(input$trait_page, {
    hits <- trait_hits()
    if (is.null(hits) || nrow(hits) == 0) return()
    npages <- max(1, ceiling(nrow(hits) / PAGE_SIZE))
    trait_page(max(1, min(as.integer(input$trait_page$n), npages)))
  }, ignoreInit = TRUE)

  observeEvent(input$trait_click, {
    gene_hits(NULL)
    trait_hits(NULL)

    # gene-card clicks carry the source gene to star; trait searches send gene:''
    gid <- input$trait_click$gene
    highlight_gene(if (is.null(gid) || !nzchar(gid)) NULL else gid)

    selected_trait(input$trait_click$trait)
    updateNavbarPage(session, "main_navbar", selected = "Results")
  })
  
  # ----------------------- code that renders the output results -------------
  
  output$search_results <- renderUI({
    mode <- input$search_type
    if (is.null(mode)) mode <- "gene"
    if (mode == "gene") {
      render_gene_results(gene_hits(), gene_page(), gstate())
    } else if (mode == "trait") {
      render_trait_results(trait_hits(), trait_page())
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
  
}

shinyApp(ui = ui, server = server)