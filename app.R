library(shiny)
library(plotly)
library(bslib)
library(DT)
library(DBI)
library(RSQLite)
source("scripts/plot_section.R")
source("scripts/search_section.R")

# Gene-based Bonferroni denominator: number of genes tested (~18,645, constant
# across traits). Read once from a representative trait file — matches the trait
# Manhattan threshold 0.05 / (genes tested).
N_GENES <- sum(!is.na(readRDS("data/Q_100001.rds")$P_mBATcombo))

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

  /* Tabs: Home, Trait, Gene, Search */
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
  
  nav_panel("Trait",
            uiOutput("trait_title"),
            card(card_header("Manhattan Plot — P_mBATcombo"), plotlyOutput("manhattan_mBATcombo")),
            card(card_header("Manhattan Plot — P_mBAT (top) vs P_fastBAT (bottom)"),
                 plotlyOutput("manhattan_mirror")),
            card(
              card_header(
                "Top Genes",
                downloadButton("download_csv", "Download CSV", class = "btn-sm btn-outline-primary float-end")
              ),
              DTOutput("gene_table")
            )
  ),
  
  nav_panel("Gene",
            div(class = "container mt-5",
                uiOutput("gene_header"),
                card(card_header("Phenome-wide associations — Phecodes (diseases)"),
                     plotlyOutput("gene_phecode_plot")),
                card(card_header("Phenome-wide associations — UK Biobank fields"),
                     plotlyOutput("gene_field_plot"))
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
                      tags$option(value = "gene",  "Gene"),
                      tags$option(value = "trait", "Trait")
                    ),
                    tags$input(
                      id          = "gene_query",
                      type        = "text",
                      class       = "form-control",
                      placeholder = "Search gene or trait..."
                    )
                ),
                # to check for the enter from the search box.
                tags$script(HTML("document.getElementById('gene_query').addEventListener('keydown', function(e)
                {if (e.key === 'Enter') {e.preventDefault();
                Shiny.setInputValue('search_submit', {query: this.value,nonce: Math.random()}, {priority: 'event'});}});")),
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
  selected_gene  <- reactiveVal("ENSG00000140718")  # default gene: FTO
  highlight_gene <- reactiveVal(NULL)   # ensembl id to star on the Trait Manhattan plots
  highlight_trait <- reactiveVal(NULL)  # trait id to star on the Gene PheWAS plots
  
  # Cross-navigation between the Trait and Gene tabs (dot clicks). Each direction
  # stars where you came from and clears the other star.
  open_trait_from_gene <- function(tid) {          # Gene PheWAS dot -> Trait tab
    if (is.null(tid) || !nzchar(tid)) return()
    highlight_gene(selected_gene()); highlight_trait(NULL)
    selected_trait(tid)
    updateNavbarPage(session, "main_navbar", selected = "Trait")
  }
  open_gene_from_trait <- function(ens) {          # Trait Manhattan dot -> Gene tab
    if (is.null(ens) || !nzchar(ens)) return()
    highlight_trait(selected_trait()); highlight_gene(NULL)
    selected_gene(ens)
    updateNavbarPage(session, "main_navbar", selected = "Gene")
  }
  
  # ---- Gene tab: two PheWAS plots for the selected gene ----
  # Fixed left-to-right order of the 18 PheWAS phecode categories (user-specified).
  PHECODE_CATS <- c("digestive", "genitourinary", "neoplasms", "musculoskeletal",
                    "circulatory system", "symptoms", "injuries & poisonings", "sense organs",
                    "respiratory", "dermatologic", "neurological", "endocrine/metabolic",
                    "infectious diseases", "mental disorders", "pregnancy complications",
                    "hematopoietic", "congenital anomalies", "other/health-services")
  
  # All of the selected gene's trait associations (ascending mBAT-combo p),
  # joined to category metadata. reactive() => computed once per gene change.
  gene_data <- reactive({
    gid <- selected_gene()
    if (is.null(gid) || !nzchar(gid)) return(NULL)
    gene_associations(gid)
  })
  
  # Shared Gene-page y-scale: max capped -log10(P_mBATcombo) over all the gene's
  # traits (both PheWAS plots), floored by the Bonferroni line, +5% headroom.
  gene_ymax <- reactive({
    df <- gene_data()
    if (is.null(df) || !nrow(df)) return(NULL)
    dmax <- max(pmin(-log10(df$P_mBATcombo), PVAL_CAP), na.rm = TRUE)
    max(dmax, -log10(0.05 / N_GENES)) * 1.05
  })
  
  output$gene_header <- renderUI({
    gid <- selected_gene()
    g   <- gene_row(gid)
    sym <- if (!is.na(g$symbol) && nzchar(g$symbol)) g$symbol else gid
    div(h2(sym),
        tags$span(gid, style = "color:#6c757d; font-size:0.9em;"),
        if (!is.na(g$genename) && nzchar(g$genename))
          tags$p(g$genename, style = "color:#6c757d; margin-top:4px;"))
  })
  
  # Phecode plot: x = the 18 PheWAS disease categories.
  output$gene_phecode_plot <- renderPlotly({
    df <- gene_data()
    if (is.null(df) || !nrow(df)) return(plotly_empty())
    d <- df[df$is_phecode == 1, , drop = FALSE]
    d$cat <- d$category; d$label <- d$category
    phewas_plot(d, levels = PHECODE_CATS, xlab = "PheWAS disease category",
                source = "gene_phewas", highlight = highlight_trait(),
                bonf = 0.05 / N_GENES, ymax = gene_ymax())
  })
  
  # Field plot: x = secondary node, grouped (and coloured) by root.
  output$gene_field_plot <- renderPlotly({
    df <- gene_data()
    if (is.null(df) || !nrow(df)) return(plotly_empty())
    d <- df[df$is_phecode == 0, , drop = FALSE]
    d$cat   <- paste(d$category_group, d$category, sep = " ▸ ")  # keep same-named secondaries distinct
    d$label <- d$category
    lv <- unique(d[order(d$category_group, d$category), c("category_group", "category")])
    levels <- paste(lv$category_group, lv$category, sep = " ▸ ")
    phewas_plot(d, levels = levels, xlab = "UK Biobank field category (grouped by root)",
                source = "gene_phewas", highlight = highlight_trait(),
                bonf = 0.05 / N_GENES, ymax = gene_ymax())
  })
  
  output$trait_title <- renderUI({
    tid   <- selected_trait()
    trait_meta <- readRDS("data/trait_metadata.rds")
    i     <- match(tid, trait_meta$trait_id)
    tname <- if (is.na(i) || is.na(trait_meta$trait_name[i])) tid else trait_meta$trait_name[i]
    add   <- if (!is.na(i) && "additional" %in% names(trait_meta)) trait_meta$additional[i] else ""
    rm(trait_meta)
    
    # `additional` is a newline-separated "Label: value" block (see build_metadata.R);
    # render each line, linking any http(s) value (download / source).
    add_ui <- NULL
    if (!is.null(add) && !is.na(add) && nzchar(add)) {
      lines <- strsplit(add, "\n", fixed = TRUE)[[1]]
      items <- lapply(lines, function(ln) {
        p   <- regexpr(": ", ln, fixed = TRUE)            # split at the FIRST ": "
        lab <- substr(ln, 1, p - 1)
        val <- substr(ln, p + 2, nchar(ln))
        vui <- if (grepl("^https?://", val)) tags$a(href = val, target = "_blank", val) else val
        tags$div(style = "padding:2px 0;", tags$strong(paste0(lab, ": ")), vui)
      })
      add_ui <- card(card_header("Additional information"), items)
    }
    
    div(
      h2(tname),
      tags$span(tid, style = "color:#6c757d; font-size:0.9em;"),
      add_ui
    )
  })
  
  # Shared Trait-page y-scale: max capped -log10(p) across all three methods,
  # floored by the Bonferroni line, +5% headroom.
  trait_ymax <- reactive({
    d <- readRDS(file.path("data", paste0(selected_trait(), ".rds")))
    dmax <- max(pmin(-log10(c(d$P_mBATcombo, d$P_mBAT, d$P_fastBAT)), PVAL_CAP), na.rm = TRUE)
    max(dmax, -log10(0.05 / N_GENES)) * 1.05
  })
  
  output$manhattan_mBATcombo <- renderPlotly({
    plot(file.path("data", paste0(selected_trait(), ".rds")), "P_mBATcombo",
         highlight_gene(), source = "manhattan", ymax = trait_ymax())
  })
  
  output$manhattan_mirror <- renderPlotly({
    mirror_plot(file.path("data", paste0(selected_trait(), ".rds")),
                "P_mBAT", "P_fastBAT", highlight_gene(), source = "manhattan",
                ymax = trait_ymax())
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
  
  # ---- dot clicks: Gene PheWAS dot -> Trait tab; Trait Manhattan dot -> Gene tab
  observeEvent(event_data("plotly_click", source = "gene_phewas"), {
    d <- event_data("plotly_click", source = "gene_phewas")
    if (!is.null(d$customdata)) open_trait_from_gene(d$customdata[[1]])
  })
  observeEvent(event_data("plotly_click", source = "manhattan"), {
    d <- event_data("plotly_click", source = "manhattan")
    if (!is.null(d$customdata)) open_gene_from_trait(d$customdata[[1]])
  })
  
  #---------------------------
  # section for search
  #---------------------------
  
  
  gene_query_d <- debounce(reactive(input$gene_query), 300)
  
  # Enter-key result lists (keyword matches held for rendering); one is non-NULL
  # at a time depending on the search mode. Each has its own page cursor.
  gene_hits  <- reactiveVal(NULL)   # data.frame(id, symbol, genename)
  trait_hits <- reactiveVal(NULL)   # data.frame(id, trait_name)
  gene_page  <- reactiveVal(1)      # current page of the gene result list
  trait_page <- reactiveVal(1)      # current page of the trait result list
  
  # Live drop-down suggestions (grep only), one per mode.
  output$gene_suggestions <- renderUI({
    if (is.null(input$search_type) || input$search_type == "trait") return(NULL)
    q <- gene_query_d()
    if (nchar(q) < 2) return(NULL)
    div(
      style = "border:1px solid #dee2e6; border-radius:4px; max-height:260px; overflow-y:auto; background:white;",
      search_results_rows(search_gene(q), "gene")
    )
  })
  
  output$trait_suggestions <- renderUI({
    if (is.null(input$search_type) || input$search_type != "trait") return(NULL)
    q <- gene_query_d()
    if (nchar(q) < 2) return(NULL)
    div(
      style = "border:1px solid #dee2e6; border-radius:4px; max-height:260px; overflow-y:auto; background:white;",
      search_results_rows(search_trait(q), "trait")
    )
  })
  
  # Enter in the search box: keyword match, held for the persistent result list.
  observeEvent(input$search_submit, {
    q <- input$search_submit$query
    if (is.null(q)) q <- ""
    q <- trimws(q)
    if (nchar(q) < 2) return()
    
    updateTextInput(session, "gene_query", value = "")
    
    mode <- if (is.null(input$search_type)) "gene" else input$search_type
    if (mode == "trait") {
      trait_hits(search_trait(q)); trait_page(1); gene_hits(NULL)
    } else {
      gene_hits(search_gene(q));   gene_page(1);  trait_hits(NULL)
    }
  }, ignoreInit = TRUE)
  
  # A gene result (dropdown or list) -> select it and jump to the Gene tab.
  observeEvent(input$gene_clicked, {
    gid <- input$gene_clicked
    if (is.null(gid) || !nzchar(gid)) return()
    updateTextInput(session, "gene_query", value = "")
    gene_hits(NULL); trait_hits(NULL)
    highlight_trait(NULL)            # arriving via Search -> no came-from star
    selected_gene(gid)
    updateNavbarPage(session, "main_navbar", selected = "Gene")
  }, ignoreInit = TRUE)
  
  # A trait result (dropdown or list) -> select it and jump to the Trait tab.
  # `gene` is the ensembl id to star (set when navigating from a plot; empty from
  # a plain search).
  observeEvent(input$trait_click, {
    updateTextInput(session, "gene_query", value = "")
    gene_hits(NULL); trait_hits(NULL)
    gid <- input$trait_click$gene
    highlight_gene(if (is.null(gid) || !nzchar(gid)) NULL else gid)
    selected_trait(input$trait_click$trait)
    updateNavbarPage(session, "main_navbar", selected = "Trait")
  })
  
  # ---- page through the gene / trait result lists (clamp to valid range) ----
  observeEvent(input$gene_page, {
    hits <- gene_hits(); if (is.null(hits) || nrow(hits) == 0) return()
    npages <- max(1, ceiling(nrow(hits) / PAGE_SIZE))
    gene_page(max(1, min(as.integer(input$gene_page$n), npages)))
  }, ignoreInit = TRUE)
  
  observeEvent(input$trait_page, {
    hits <- trait_hits(); if (is.null(hits) || nrow(hits) == 0) return()
    npages <- max(1, ceiling(nrow(hits) / PAGE_SIZE))
    trait_page(max(1, min(as.integer(input$trait_page$n), npages)))
  }, ignoreInit = TRUE)
  
  # ----------------------- render the Enter-key result list -------------
  output$search_results <- renderUI({
    mode <- input$search_type
    if (is.null(mode)) mode <- "gene"
    if (mode == "gene") render_results(gene_hits(),  "gene",  gene_page())
    else                render_results(trait_hits(), "trait", trait_page())
  })
  
  # ----------------
  # home page section
  #-----------------
  
  observeEvent(input$home_browse, {
    updateNavbarPage(session, "main_navbar", selected = "Trait")
  })
  observeEvent(input$home_search, {
    updateNavbarPage(session, "main_navbar", selected = "Search")
  })
  
}

shinyApp(ui = ui, server = server)