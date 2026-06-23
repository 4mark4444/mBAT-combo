#----------------------SEARCH SECTION. ----------------------

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


#---------------------- server helper functions -----------------------------

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



# This is the function that randers the UI for the drop down suggestions. 
search_results_rows <- function(search_results, type){
  rows <- list()
  for (i in seq_len(nrow(search_results))){
    if (type == "gene"){
      event <- "gene_clicked"
      title <- search_results$symbol[i]
    }else{
      event <- "trait_click"
      title <- search_results$trait_name[i]
    }
    rows[[i]] <- div(
      style   = "padding:8px 12px; border-bottom:1px solid #f0f0f0; cursor:pointer;",
      onclick = sprintf("Shiny.setInputValue('%s', '%s', {priority:'event'})", event, search_results$id[i]),
      tags$b(title),
      tags$span(search_results$id[i], style = "color:#6c757d; font-size:0.85em; margin-left:6px;"),
      if (type == "gene") tags$br(),
      if (type == "gene") tags$span(search_results$genename[i], style = "color:#6c757d; font-size:0.8em;")
    )
  }
  
  return (rows)
}


# The inverted index search contains ONLY the id. Hence to make the gene actually human reable, we need to fetch the name. 
fetch_names <- function (results){
  gene_map <- readRDS("data/gene_symbol_map.rds")
  trait_meta <- readRDS("data/trait_metadata.rds")
  
  results$symbol     <- gene_map$SYMBOL[match(results$ensembl_id, gene_map$ENSEMBL)]
  results$trait_name <- trait_meta$trait_name[match(results$trait_id, trait_meta$trait_id)]
  
  rm(gene_map)
  rm(trait_meta)
  return (results)
}


# Given a gene this function creates the UI that add the asociated traits to it. 
association_card <- function(g) {
  ensembl_id <- g$ensembl_id[1]
  symbol     <- if (!is.na(g$symbol[1])) g$symbol[1] else ensembl_id
  
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
  
  div(
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

# This is the UI render function for the search results. After the user press enter, this is the results that they see. 
render_results <- function(data, mode) {
  if (is.null(data) || nrow(data) == 0) return(NULL)
  
  cards <- list()
  
  if (mode == "gene") {
    data        <- fetch_names(data)
    gene_groups <- split(data, data$ensembl_id)
    for (id in names(gene_groups)) {
      cards[[id]] <- association_card(gene_groups[[id]])
    }
  } else if (mode == "trait") {
    for (i in seq_len(nrow(data))) {
      cards[[i]] <- div(
        class = "card mb-2",
        div(class = "card-body py-2",
            style   = "cursor:pointer;",
            onclick = sprintf("Shiny.setInputValue('trait_click', '%s', {priority:'event'})", data$id[i]),
            tags$b(data$trait_name[i]),
            tags$span(data$id[i], style = "color:#6c757d; font-size:0.85em; margin-left:8px;")
        )
      )
    }
  }
  
  div(class = "mt-3", unname(cards))
}
