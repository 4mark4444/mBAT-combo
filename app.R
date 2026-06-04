library(shiny)
library(plotly)
library(bslib)
library(DBI)
library(RSQLite)
library(DT)

# --------------------------------------
# helper function for plotting the graph. 
# --------------------------------------

plot <- function(path, p){
  data <- readRDS(path)
  data$y <- -log10(data[[p]])
  chromo <- sort(unique(data$Chr))
  
  colors <- list(
    "P_mBATcombo" = c("#85b3d1", "#2c6fad"),
    "P_mBAT"      = c("#ff9999", "#cc0000"),
    "P_fastBAT"   = c("#99d899", "#228b22")
  )
  
  # assign color 
  color <- c()
  for (i in 1:length(chromo)){
    if (i %%2 == 0){
      color[chromo[i]] <-  colors[[p]][1] # R list does not start from 0 ?!
    }
    else{
      color[chromo[i]] <-  colors[[p]][2]
    }
  }
  data$color <- color[data$Chr]
  
  # calculate the middle point so that we can put the label here. 
  chr_mid <- c()
  for (i in chromo) {
    chr_mid[as.character(i)] <- median(data$cum_pos[data$Chr == i])
  }
  # add the hover text display. 
  data$hover <- paste0(                                                                                                                                                                                                                                                     
    "<b>", data$Symbol, "</b><br>",                                                                                                                                                                                                                                       
    "Ensembl: ",        data$EnsemblID, "<br>",                                                                                                                                                                                                                             
    "Chr",              data$Chr, ":", format(data$Start, big.mark = ","),                                                                                                                                                                                                   
    "–",                format(data$End, big.mark = ","), "<br>",                                                                                                                                                                                                           
    "No.SNPs: ",        data$No.SNPs, "<br>",                                                                                                                                                                                                                               
    "P_mBATcombo: ",    formatC(data$P_mBATcombo,   format = "e", digits = 3), "<br>",                                                                                                                                                                                      
    "P_mBAT: ",         formatC(data$P_mBAT,        format = "e", digits = 3), "<br>",                                                                                                                                                                                      
    "P_fastBAT: ",      formatC(data$P_fastBAT,     format = "e", digits = 3), "<br>",                                                                                                                                                                                      
    "TopSNP: ",         data$TopSNP, "<br>",                                                                                                                                                                                                                                
    "TopSNP p: ",       formatC(data$TopSNP_Pvalue, format = "e", digits = 3), "<br>",                                                                                                                                                                                      
    "No.Eigenvalues: ", data$No.Eigenvalues                                                                                                                                                                                                                                 
  )
  
  graph <- plot_ly(data = data, x = ~cum_pos, y = ~y, 
                   type = "scatter", mode = "markers", marker = list(color = ~color), 
                   text = ~hover, hoverinfo = "text"
                   )
  graph <- layout(graph, 
                  xaxis = list(title = "Chromosome", tickvals = chr_mid, ticktext = names(chr_mid)),  
                  yaxis = list(title = paste0("-log10(", p, ")"))
                  )
  rm(data)                
  return (graph)
}

# ---------------------------
# helper function for the search. This is the stage 1 of the search function. 
# ---------------------------
search_gene <- function(query) {
  gene_map <- readRDS("data/gene_symbol_map.rds")
  hits <- gene_map[
    grepl(query, gene_map$SYMBOL,   ignore.case = TRUE) |
    grepl(query, gene_map$GENENAME, ignore.case = TRUE) |
    grepl(query, gene_map$ENSEMBL, ignore.case = TRUE), 
  ]
  rm(gene_map)
  return( data.frame(
    id       = hits$ENSEMBL,
    symbol   = hits$SYMBOL,
    genename = hits$GENENAME,
    stringsAsFactors = FALSE
  ))
}

# ------------------------------------------
# Inverted index search using SQL. 
# ------------------------------------------

reverse_index <- function(id_list, p1, p2, p3) {
  if (length(id_list) == 0) return(data.frame())
  
  con <- dbConnect(RSQLite::SQLite(), "data/search.sqlite")
  on.exit(dbDisconnect(con))
  
  placeholders <- paste(rep("?", length(id_list)), collapse = ", ")
  sql <- sprintf(
    "SELECT ensembl_id, trait_id, P_mBATcombo, P_mBAT, P_fastBAT
     FROM associations
     WHERE ensembl_id IN (%s)
       AND P_mBATcombo < ?
       AND P_mBAT      < ?
       AND P_fastBAT   < ?",
    placeholders
  )
  
  dbGetQuery(con, sql, params = c(as.list(id_list), p1, p2, p3))
}

#-------------------------
# helper function for search trait. 
#-------------------------

search_trait <- function(query) {
  trait_meta <- readRDS("data/trait_metadata.rds")
  hits <- trait_meta[
    grepl(query, trait_meta$trait_name, ignore.case = TRUE),
  ]
  rm(trait_meta)
  return (data.frame(
    id         = hits$trait_id,
    trait_name = hits$trait_name,
    stringsAsFactors = FALSE
  ))
}

#----------------------------
# Helper function for webUI. 
#----------------------------

# Drop down for the filter
filter_row <- function(label, id_cb, id_num) {
  return ( 
    div(
    style = "display:flex; align-items:center; gap:10px; margin-bottom:8px;",
    tags$input(type = "checkbox", id = id_cb, checked = NA,
               style = "width:16px; height:16px; cursor:pointer;"),
    tags$span(label, style = "width:140px;"),
    tags$span("<"),
    tags$input(type = "number", id = id_num, value = "0.05",
               min = "0", max = "1", step = "0.001",
               autocomplete = "off",
               class = "form-control form-control-sm",
               style = "width:90px;")
  ) 
  )
}

# drop down search results

search_result_gene_rows <- function(search_result){
  rows <- list()
  for (i in seq_len(nrow(search_result))) {
    gene_id <- search_result$id[i]
    rows[[i]] <- div(
      style = "padding:8px 12px; border-bottom:1px solid #f0f0f0; cursor:pointer;", 
      onclick = sprintf("Shiny.setInputValue('gene_clicked', '%s', {priority:'event'})", gene_id),
      tags$b(search_result$symbol[i]),
      tags$span(
        search_result$id[i],
        style = "color:#6c757d; font-size:0.85em; margin-left:6px;"
      ),
      tags$br(),
      tags$span(
        search_result$genename[i],
        style = "color:#6c757d; font-size:0.8em;"
      )
    )
  }
  
  return (rows)
}

# Search trait drop down. TODO: merge the two drop down function. 

search_result_trait_rows <- function(search_result){
  rows <- list()
  for (i in seq_len(nrow(search_result))) {
    trait_id <- search_result$id[i]
    rows[[i]] <- div(
      style = "padding:8px 12px; border-bottom:1px solid #f0f0f0; cursor:pointer;",
      onclick = sprintf("Shiny.setInputValue('trait_click', '%s', {priority:'event'})", trait_id),
      tags$b(search_result$trait_name[i]),
      tags$span(
        search_result$id[i],
        style = "color:#6c757d; font-size:0.85em; margin-left:6px;"
      )
    )
  }
  
  return (rows)
}


stat_box <- function(label, value) {
  div(
    style = "color: white; text-align: center;",
    div(style = "font-size: 1.6rem; font-weight: 700;", value),
    div(style = "font-size: 0.8rem; opacity: 0.75; text-transform: uppercase; letter-spacing: 0.05em;", label)
  )
}

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
      search_result_gene_rows(search_gene(q))
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
      search_result_trait_rows(search_trait(q))
    )
    
  })
  
  search_results_data <- reactiveVal(NULL)
  trait_search_data <- reactiveVal(NULL) 
  
  observeEvent(input$gene_clicked, {
    
    gene_id <- input$gene_clicked
    if (is.null(gene_id) || !nzchar(gene_id)) return()
    
    updateTextInput(session, "gene_query", value = "")
    
    p1 <- if (isTRUE(input$cb_p1)) input$thresh_p1 else 1
    p2 <- if (isTRUE(input$cb_p2)) input$thresh_p2 else 1
    p3 <- if (isTRUE(input$cb_p3)) input$thresh_p3 else 1
    
    results <- reverse_index(gene_id, p1, p2, p3)
    search_results_data(results)
  }, ignoreInit = TRUE)
  
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
      search_results_data(NULL)
      
      return()
    }
    
    genes <- search_gene(q)
    if (nrow(genes) == 0) {
      search_results_data(data.frame())
      trait_search_data(NULL)
      return()
    }
    
    p1 <- if (isTRUE(input$cb_p1)) input$thresh_p1 else 1
    p2 <- if (isTRUE(input$cb_p2)) input$thresh_p2 else 1
    p3 <- if (isTRUE(input$cb_p3)) input$thresh_p3 else 1
    
    results <- reverse_index(genes$id, p1, p2, p3)
    
    search_results_data(results)
    trait_search_data(NULL)
  }, ignoreInit = TRUE)
  
  output$search_results <- renderUI({
    
    traits <- trait_search_data()
    if (!is.null(traits)) {
      if (nrow(traits) == 0) return(NULL)
      
      return(
        div(
          class = "mt-3",
          style = "border:1px solid #dee2e6; border-radius:4px; max-height:260px; overflow-y:auto; background:white;",
          search_result_trait_rows(traits)
        )
      )
    }
    
    results <- search_results_data()
    if (is.null(results) || nrow(results) == 0) return(NULL)
    
    gene_map <- readRDS("data/gene_symbol_map.rds")
    trait_meta <- readRDS("data/trait_metadata.rds")
    
    results$symbol     <- gene_map$SYMBOL[match(results$ensembl_id, gene_map$ENSEMBL)]
    results$trait_name <- trait_meta$trait_name[match(results$trait_id, trait_meta$trait_id)]
    
    rm(gene_map)
    rm(trait_meta)
    
    gene_groups <- split(results, results$ensembl_id)
    
    cards <- list()
    for (ensembl_id in names(gene_groups)) {
      g      <- gene_groups[[ensembl_id]]
      symbol <- if (!is.na(g$symbol[1])) g$symbol[1] else ensembl_id
      
      trait_rows <- list()
      for (i in seq_len(nrow(g))) {
        tname <- if (!is.na(g$trait_name[i])) g$trait_name[i] else g$trait_id[i]
        trait_rows[[i]] <- tags$li(
          style = "padding:4px 0;",
          tags$a(
            href    = "#",
            onclick = sprintf("Shiny.setInputValue('trait_click', '%s', {priority:'event'})", g$trait_id[i]),
            tname
          ),
          tags$small(
            style = "color:#6c757d; margin-left:16px;",
            sprintf("P_mBATcombo: %s  |  P_mBAT: %s  |  P_fastBAT: %s",
                    formatC(g$P_mBATcombo[i], format = "e", digits = 2),
                    formatC(g$P_mBAT[i],      format = "e", digits = 2),
                    formatC(g$P_fastBAT[i],   format = "e", digits = 2))
          )
        )
      }
      
      cards[[ensembl_id]] <- div(
        class = "card mb-3",
        div(class = "card-header",
            tags$b(symbol),
            tags$span(ensembl_id, style = "color:#6c757d; font-size:0.85em; margin-left:8px;")
        ),
        div(class = "card-body py-2",
            tags$ul(style = "margin:0; padding-left:20px;", trait_rows)
        )
      )
    }
    
    div(class = "mt-3", cards)
  })
  
  observeEvent(input$trait_click, {
    trait_search_data(NULL)
    search_results_data(NULL)
    
    selected_trait(input$trait_click)
    updateNavbarPage(session, "main_navbar", selected = "Results")
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