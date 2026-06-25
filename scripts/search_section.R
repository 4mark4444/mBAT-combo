#----------------------SEARCH SECTION. ----------------------

# Page size shared by gene-card pagination, per-card "load more", and trait-card
# pagination. One screenful is 10 of everything.
PAGE_SIZE <- 10

# ---------------------------
# Stage 1 text match (in-memory grep over the RDS maps). Powers both the live
# dropdown suggestions and the Enter-key submit. No SQLite here.
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
# clicked, so the gene-card path and the Enter path share the same shape).
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
# Stage 2: lazy, per-gene, rowid-cursor read of the inverted index.
# ------------------------------------------
# The rebuilt search.sqlite stores rows physically ordered (ensembl_id,
# P_mBATcombo), so within a gene rowid order IS ascending mBAT-combo p-value. The
# index on ensembl_id is really (ensembl_id, rowid), so "ensembl_id = ? AND rowid >
# ?" is a direct seek -> no run-time sort, no re-walking already-shown rows. We ask
# for one extra row (LIMIT n+1) so the caller can tell whether a "Load more" remains.
#
# p_* are p-value thresholds; NA means that filter is inactive (so we don't emit the
# clause at all -- keeps the default "no filter" case index-only and fast).

.pfilter <- function(p_combo, p_mbat, p_fastbat) {
  cl <- character(0); pr <- list()
  if (!is.null(p_combo)   && !is.na(p_combo))   { cl <- c(cl, "P_mBATcombo < ?"); pr <- c(pr, p_combo) }
  if (!is.null(p_mbat)    && !is.na(p_mbat))    { cl <- c(cl, "P_mBAT < ?");      pr <- c(pr, p_mbat) }
  if (!is.null(p_fastbat) && !is.na(p_fastbat)) { cl <- c(cl, "P_fastBAT < ?");   pr <- c(pr, p_fastbat) }
  list(clause = if (length(cl)) paste("AND", paste(cl, collapse = " AND ")) else "",
       params = pr)
}

reverse_index_page <- function(gene_id, p_combo, p_mbat, p_fastbat,
                               after_rowid = 0, limit = PAGE_SIZE) {
  con <- dbConnect(RSQLite::SQLite(), "data/search.sqlite")
  on.exit(dbDisconnect(con))

  f <- .pfilter(p_combo, p_mbat, p_fastbat)
  sql <- sprintf(
    "SELECT rowid AS rid, trait_id, P_mBATcombo, P_mBAT, P_fastBAT
       FROM associations
      WHERE ensembl_id = ? AND rowid > ? %s
      ORDER BY rowid
      LIMIT ?",
    f$clause
  )
  dbGetQuery(con, sql, params = c(list(gene_id, after_rowid), f$params, list(limit)))
}


#---------------------- UI render helpers -----------------------------

# Drop-down filter row (one p-value). Default UNCHECKED == filter inactive, so a
# plain search returns everything (the user's "keep current behavior").
filter_row <- function(label, id_cb, id_num) {
  return (
    div(
      style = "display:flex; align-items:center; gap:10px; margin-bottom:8px;",
      tags$input(type = "checkbox", id = id_cb,
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

# Rows for the live drop-down suggestions (Stage 1 only).
search_results_rows <- function(search_results, type){
  rows <- list()
  for (i in seq_len(nrow(search_results))){
    if (type == "gene"){
      title   <- search_results$symbol[i]
      onclick <- sprintf("Shiny.setInputValue('gene_clicked', '%s', {priority:'event'})",
                         search_results$id[i])
    }else{
      title   <- search_results$trait_name[i]
      # trait dropdown carries no gene to highlight
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

# One gene card: header + the traits loaded so far + (optional) "Load more".
# `rows` must already carry a trait_name column.
gene_card <- function(ensembl_id, symbol, rows, more) {
  symbol <- if (!is.na(symbol) && nzchar(symbol)) symbol else ensembl_id

  if (is.null(rows) || nrow(rows) == 0) {
    body <- tags$div(class = "text-muted", "No associations pass the current filters.")
  } else {
    items <- lapply(seq_len(nrow(rows)), function(i) {
      tname <- if (!is.na(rows$trait_name[i])) rows$trait_name[i] else rows$trait_id[i]
      tags$li(
        style = "padding:4px 0;",
        tags$a(
          href    = "#",
          onclick = sprintf(
            "Shiny.setInputValue('trait_click', {trait:'%s', gene:'%s'}, {priority:'event'})",
            rows$trait_id[i], ensembl_id),
          tname
        ),
        tags$small(
          style = "color:#6c757d; margin-left:16px;",
          sprintf("P_mBATcombo: %s  |  P_mBAT: %s  |  P_fastBAT: %s",
                  formatC(rows$P_mBATcombo[i], format = "e", digits = 2),
                  formatC(rows$P_mBAT[i],      format = "e", digits = 2),
                  formatC(rows$P_fastBAT[i],   format = "e", digits = 2))
        )
      )
    })
    body <- tags$ul(style = "margin:0; padding-left:20px;", items)
  }

  load_more <- if (isTRUE(more)) tags$button(
    class   = "btn btn-sm btn-outline-primary mt-2",
    onclick = sprintf("Shiny.setInputValue('load_more', {id:'%s', nonce:Math.random()}, {priority:'event'})", ensembl_id),
    "Load more"
  ) else NULL

  div(
    class = "card mb-3",
    div(class = "card-header",
        tags$b(symbol),
        tags$span(ensembl_id, style = "color:#6c757d; font-size:0.85em; margin-left:8px;"),
        if (!is.null(rows) && nrow(rows) > 0)
          tags$span(sprintf("%d shown", nrow(rows)),
                    style = "color:#6c757d; font-size:0.8em; margin-left:8px;")
    ),
    div(class = "card-body py-2", body, load_more)
  )
}

# Gene-search results: one page (<=PAGE_SIZE) of gene cards + the page bar.
# `state` is the named list (one entry per loaded gene): list(rows, cursor, more).
render_gene_results <- function(hits, page, state) {
  if (is.null(hits)) return(NULL)
  if (nrow(hits) == 0) return(div(class = "mt-3 text-muted", "No matching genes."))

  trait_meta <- readRDS("data/trait_metadata.rds")
  npages <- max(1, ceiling(nrow(hits) / PAGE_SIZE))
  page   <- max(1, min(page, npages))
  idx    <- ((page - 1) * PAGE_SIZE + 1):min(page * PAGE_SIZE, nrow(hits))

  cards <- lapply(idx, function(i) {
    id   <- hits$id[i]
    st   <- state[[id]]
    rows <- if (is.null(st)) NULL else st$rows
    if (!is.null(rows) && nrow(rows) > 0) {
      rows$trait_name <- trait_meta$trait_name[match(rows$trait_id, trait_meta$trait_id)]
    }
    gene_card(id, hits$symbol[i], rows, more = !is.null(st) && isTRUE(st$more))
  })
  rm(trait_meta)

  div(class = "mt-3", cards, page_bar(page, npages, "gene_page"))
}

# Trait-search results: one page of clickable trait cards + the page bar. No DB.
render_trait_results <- function(hits, page) {
  if (is.null(hits)) return(NULL)
  if (nrow(hits) == 0) return(div(class = "mt-3 text-muted", "No matching traits."))

  npages <- max(1, ceiling(nrow(hits) / PAGE_SIZE))
  page   <- max(1, min(page, npages))
  idx    <- ((page - 1) * PAGE_SIZE + 1):min(page * PAGE_SIZE, nrow(hits))

  cards <- lapply(idx, function(i) {
    div(
      class = "card mb-2",
      div(class = "card-body py-2",
          style   = "cursor:pointer;",
          onclick = sprintf("Shiny.setInputValue('trait_click', {trait:'%s', gene:''}, {priority:'event'})", hits$id[i]),
          tags$b(hits$trait_name[i]),
          tags$span(hits$id[i], style = "color:#6c757d; font-size:0.85em; margin-left:8px;")
      )
    )
  })

  div(class = "mt-3", cards, page_bar(page, npages, "trait_page"))
}
