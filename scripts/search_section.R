#----------------------SEARCH SECTION. ----------------------
#
# Plain keyword search only. Both gene and trait search are in-memory greps over
# the RDS maps; there is no SQLite / inverted-index / p-value filtering here any
# more. A search result's only job is to navigate: gene rows fire `gene_clicked`
# (-> Gene tab), trait rows fire `trait_click` (-> Trait tab).

# Enter-key result list page size (one screenful).
PAGE_SIZE <- 10

# ---------------------------
# Text match (in-memory grep over the RDS maps). Powers both the live dropdown
# suggestions and the Enter-key submit list.
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

# Resolve one ensembl id -> a one-row hits frame (used when a gene suggestion is
# clicked, so the click path and the Enter path share the same shape).
gene_row <- function(ensembl_id) {
  gene_map <- readRDS("data/gene_symbol_map.rds")
  i <- match(ensembl_id, gene_map$ENSEMBL)
  out <- data.frame(
    id       = ensembl_id,
    symbol   = if (is.na(i)) NA_character_ else gene_map$SYMBOL[i],
    genename = if (is.na(i)) NA_character_ else gene_map$GENENAME[i],
    stringsAsFactors = FALSE
  )
  rm(gene_map)
  out
}

# ------------------------------------------
# Gene-centric read of the inverted index (drives the Gene tab's PheWAS plots).
# ------------------------------------------
# One indexed seek returns ALL of a gene's trait associations. search.sqlite is
# physically ordered (ensembl_id, P_mBATcombo), so rowid order == ascending
# mBAT-combo p -> the rows come back most-significant-first with no run-time sort,
# which is exactly the within-category ranking the plots want. Each trait is then
# joined (vectorised match(), no loop) to its category / phecode flag in
# trait_metadata.rds.
gene_associations <- function(gene_id) {
  con <- dbConnect(RSQLite::SQLite(), "data/search.sqlite")
  on.exit(dbDisconnect(con))
  df <- dbGetQuery(con,
                   "SELECT trait_id, P_mBATcombo, P_mBAT, P_fastBAT
       FROM associations WHERE ensembl_id = ? ORDER BY rowid",
                   params = list(gene_id))
  if (!nrow(df)) return(df)
  meta <- readRDS("data/trait_metadata.rds")
  i <- match(df$trait_id, meta$trait_id)
  df$trait_name     <- meta$trait_name[i]
  df$is_phecode     <- meta$is_phecode[i]
  df$category       <- meta$category[i]
  df$category_group <- meta$category_group[i]
  rm(meta)
  df
}

#---------------------- UI render helpers -----------------------------

# Clickable rows for both the live drop-down suggestions and the Enter-key result
# list. Gene rows navigate to the Gene tab (`gene_clicked`); trait rows navigate
# to the Trait tab (`trait_click`, no gene to highlight from search).
search_results_rows <- function(search_results, type){
  rows <- list()
  for (i in seq_len(nrow(search_results))){
    if (type == "gene"){
      title   <- search_results$symbol[i]
      onclick <- sprintf("Shiny.setInputValue('gene_clicked', '%s', {priority:'event'})",
                         search_results$id[i])
    }else{
      title   <- search_results$trait_name[i]
      onclick <- sprintf("Shiny.setInputValue('trait_click', {trait:'%s', gene:''}, {priority:'event'})",
                         search_results$id[i])
    }
    rows[[i]] <- div(
      style   = "padding:8px 12px; border-bottom:1px solid #f0f0f0; cursor:pointer;",
      onclick = onclick,
      tags$b(title),
      tags$span(search_results$id[i], style = "color:#6c757d; font-size:0.85em; margin-left:6px;"),
      if (type == "gene") tags$br(),
      if (type == "gene") tags$span(search_results$genename[i], style = "color:#6c757d; font-size:0.8em;")
    )
  }
  return (rows)
}

# Prev / "Page X of Y" / Next bar. Emits Shiny.setInputValue(<event>, {n, nonce}).
page_bar <- function(page, npages, event) {
  if (npages <= 1) return(NULL)
  nav_btn <- function(label, target, disabled) {
    tags$button(
      class = sprintf("btn btn-sm btn-outline-secondary%s", if (disabled) " disabled" else ""),
      onclick = if (disabled) NULL else
        sprintf("Shiny.setInputValue('%s', {n:%d, nonce:Math.random()}, {priority:'event'})", event, target),
      label
    )
  }
  div(
    style = "display:flex; align-items:center; gap:12px; margin-top:8px;",
    nav_btn("Prev", page - 1, page <= 1),
    tags$span(sprintf("Page %d of %d", page, npages), style = "color:#6c757d;"),
    nav_btn("Next", page + 1, page >= npages)
  )
}

# Paginated clickable result list for the Enter-key search, shared by both modes.
# Slices `hits` to the current page (PAGE_SIZE rows) and appends the page bar; the
# bar's event id (`gene_page` / `trait_page`) is chosen by `mode`.
render_results <- function(hits, mode, page) {
  if (is.null(hits)) return(NULL)
  if (nrow(hits) == 0) return(div(class = "mt-3 text-muted", "No matches."))
  npages <- max(1, ceiling(nrow(hits) / PAGE_SIZE))
  page   <- max(1, min(page, npages))
  idx    <- ((page - 1) * PAGE_SIZE + 1):min(page * PAGE_SIZE, nrow(hits))
  event  <- if (mode == "gene") "gene_page" else "trait_page"
  div(class = "mt-3",
      div(style = "border:1px solid #dee2e6; border-radius:4px; background:white;",
          search_results_rows(hits[idx, , drop = FALSE], mode)),
      page_bar(page, npages, event))
}
