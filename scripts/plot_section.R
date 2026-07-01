# The helper function for plotting the plot.

PVAL_CAP <- 300   # cap -log10(p): underflow (p==0 -> Inf) renders here; extremes flatten

plot <- function(path, p, highlight = NULL, source = NULL, ymax = NULL){
  data <- readRDS(path)
  data$y <- pmin(-log10(data[[p]]), PVAL_CAP)
  
  # Bonferroni significance threshold: 0.05 / (genes tested by this method)
  n_tests <- sum(!is.na(data[[p]]))
  bonf_p  <- 0.05 / n_tests
  bonf_y  <- -log10(bonf_p)
  
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
  
  if (is.null(ymax)) ymax <- max(c(data$y, bonf_y), na.rm = TRUE) * 1.05
  
  graph <- plot_ly(data = data, x = ~cum_pos, y = ~y, source = source,
                   customdata = ~EnsemblID,
                   type = "scatter", mode = "markers", marker = list(color = ~color),
                   text = ~hover, hoverinfo = "text"
  )
  graph <- layout(graph,
                  xaxis = list(title = "Chromosome", tickvals = chr_mid, ticktext = names(chr_mid)),
                  yaxis = list(title = paste0("-log10(", p, ")"), range = c(0, ymax)),
                  shapes = if (n_tests > 0) list(list(
                    type = "line",
                    xref = "paper", x0 = 0, x1 = 1,
                    yref = "y",     y0 = bonf_y, y1 = bonf_y,
                    line = list(color = "red", dash = "dash", width = 1)
                  )),
                  annotations = if (n_tests > 0) list(list(
                    xref = "paper", x = 1, xanchor = "right",
                    yref = "y",     y = bonf_y, yanchor = "bottom",
                    text = paste0("Bonferroni p = ", formatC(bonf_p, format = "e", digits = 2)),
                    showarrow = FALSE,
                    font = list(color = "red", size = 11)
                  ))
  )
  
  # Highlight the gene the user clicked through from gene search (gold star +
  # symbol label). Skip silently if the gene has no plottable point for this method.
  if (!is.null(highlight) && nzchar(highlight)) {
    hi <- which(data$EnsemblID == highlight)
    if (length(hi) == 1 && is.finite(data$y[hi])) {
      graph <- add_trace(graph,
                         x = data$cum_pos[hi], y = data$y[hi],
                         customdata = data$EnsemblID[hi],
                         type = "scatter", mode = "markers+text", inherit = FALSE,
                         marker = list(color = "gold", size = 14, symbol = "star",
                                       line = list(color = "black", width = 1.5)),
                         text = data$Symbol[hi], textposition = "top center",
                         textfont = list(color = "black", size = 12),
                         hovertext = data$hover[hi], hoverinfo = "text",
                         showlegend = FALSE)
    }
  }
  
  graph <- event_register(graph, "plotly_click")
  rm(data)
  return (graph)
}

# ---------------------------------------------------------------------------
# Mirror ("Miami") plot: two methods on one shared x-axis. `p_top` points up
# (y = -log10 p); `p_bottom` points down (y = -(-log10 p)). y-tick labels show
# |value| so both halves read as magnitudes. Same per-chromosome colouring,
# hover, Bonferroni line, PVAL_CAP cap, gold-star highlight and click wiring as
# plot(). Colours: p_top uses its palette (mBAT reds), p_bottom its own (fastBAT
# greens), still alternating by chromosome.
# ---------------------------------------------------------------------------
mirror_plot <- function(path, p_top, p_bottom, highlight = NULL, source = NULL, ymax = NULL) {
  data  <- readRDS(path)
  y_top <-  pmin(-log10(data[[p_top]]),    PVAL_CAP)   # up
  y_bot <- -pmin(-log10(data[[p_bottom]]), PVAL_CAP)   # down (negated)
  
  chromo  <- sort(unique(data$Chr))
  chr_mid <- vapply(chromo, function(i) median(data$cum_pos[data$Chr == i]), numeric(1))
  names(chr_mid) <- as.character(chromo)
  
  colors <- list(P_mBATcombo = c("#85b3d1","#2c6fad"),
                 P_mBAT       = c("#ff9999","#cc0000"),
                 P_fastBAT    = c("#99d899","#228b22"))
  chr_col <- function(p) {
    cc <- ifelse(seq_along(chromo) %% 2 == 0, colors[[p]][1], colors[[p]][2])
    setNames(cc, as.character(chromo))[as.character(data$Chr)]
  }
  
  data$hover <- paste0(
    "<b>", data$Symbol, "</b><br>",
    "Ensembl: ",        data$EnsemblID, "<br>",
    "Chr", data$Chr, ":", format(data$Start, big.mark = ","),
    "–", format(data$End, big.mark = ","), "<br>",
    "No.SNPs: ",        data$No.SNPs, "<br>",
    "P_mBATcombo: ",    formatC(data$P_mBATcombo,   format = "e", digits = 3), "<br>",
    "P_mBAT: ",         formatC(data$P_mBAT,        format = "e", digits = 3), "<br>",
    "P_fastBAT: ",      formatC(data$P_fastBAT,     format = "e", digits = 3), "<br>",
    "TopSNP: ",         data$TopSNP, "<br>",
    "TopSNP p: ",       formatC(data$TopSNP_Pvalue, format = "e", digits = 3), "<br>",
    "No.Eigenvalues: ", data$No.Eigenvalues)
  
  # Bonferroni line per method (top positive, bottom negated)
  bt_y <-  -log10(0.05 / sum(!is.na(data[[p_top]])))
  bb_y <- -(-log10(0.05 / sum(!is.na(data[[p_bottom]]))))
  
  # shared data-driven symmetric scale + |value| tick labels
  if (is.null(ymax)) ymax <- max(abs(c(y_top, y_bot, bt_y, bb_y)), na.rm = TRUE) * 1.05
  ticks <- pretty(c(0, ymax)); ticks <- ticks[ticks <= ymax]
  tv    <- sort(unique(c(-ticks, ticks)))
  
  g <- plot_ly(source = source)
  g <- add_trace(g, x = data$cum_pos, y = y_top, customdata = data$EnsemblID,
                 type = "scatter", mode = "markers",
                 marker = list(color = chr_col(p_top)),
                 text = data$hover, hoverinfo = "text", showlegend = FALSE)
  g <- add_trace(g, x = data$cum_pos, y = y_bot, customdata = data$EnsemblID,
                 type = "scatter", mode = "markers",
                 marker = list(color = chr_col(p_bottom)),
                 text = data$hover, hoverinfo = "text", showlegend = FALSE)
  g <- layout(g,
              xaxis = list(title = "Chromosome", tickvals = chr_mid, ticktext = names(chr_mid)),
              yaxis = list(title = "-log10(p)", range = c(-ymax, ymax),
                           tickvals = tv, ticktext = as.character(abs(tv)),
                           zeroline = TRUE, zerolinecolor = "#888", zerolinewidth = 1),
              shapes = list(
                list(type="line", xref="paper", x0=0, x1=1, yref="y", y0=bt_y, y1=bt_y,
                     line=list(color="red", dash="dash", width=1)),
                list(type="line", xref="paper", x0=0, x1=1, yref="y", y0=bb_y, y1=bb_y,
                     line=list(color="red", dash="dash", width=1))),
              annotations = list(
                list(xref="paper", x=0.005, xanchor="left", yref="paper", y=0.99, yanchor="top",
                     text=paste0(p_top, " ↑"),    showarrow=FALSE, font=list(size=12, color="#333")),
                list(xref="paper", x=0.005, xanchor="left", yref="paper", y=0.01, yanchor="bottom",
                     text=paste0(p_bottom, " ↓"), showarrow=FALSE, font=list(size=12, color="#333"))))
  
  # gold-star the came-from gene on BOTH halves (symbol label on the top one)
  if (!is.null(highlight) && nzchar(highlight)) {
    hi <- which(data$EnsemblID == highlight)
    if (length(hi) == 1) {
      star <- function(g, yv, txt) add_trace(g, x = data$cum_pos[hi], y = yv,
                                             customdata = data$EnsemblID[hi], type = "scatter",
                                             mode = "markers+text", inherit = FALSE,   # empty txt renders nothing
                                             marker = list(color="gold", size=14, symbol="star", line=list(color="black", width=1.5)),
                                             text = txt, textposition = "top center", textfont = list(color="black", size=12),
                                             hovertext = data$hover[hi], hoverinfo = "text", showlegend = FALSE)
      if (is.finite(y_top[hi])) g <- star(g, y_top[hi], data$Symbol[hi])
      if (is.finite(y_bot[hi])) g <- star(g, y_bot[hi], "")
    }
  }
  
  g <- event_register(g, "plotly_click")
  rm(data)
  g
}

# ---------------------------------------------------------------------------
# Gene PheWAS plot (Gene tab). One marker per trait, laid out along x in bands
# by category; y = -log10(P_mBATcombo). `df` is expected pre-sorted ascending
# P_mBATcombo (search.sqlite's physical order) so the stable sort below keeps
# each category band most-significant-first "for free".
#
# `df` must carry: trait_name, P_mBATcombo, P_mBAT, P_fastBAT, and
#   cat   = x grouping key (its factor levels come from `levels`)
#   label = x-tick text to show for its category
# `levels` fixes the left-to-right category order. Each category gets an equal-
# width slot (so tick labels are evenly spaced and never collide, regardless of
# how many traits fall in it) and its own rainbow colour; within a slot the
# points are spread in rank order (most-significant on the left).
#
# `source` wires plotly click events (each point's customdata = its trait_id);
# `highlight` is a trait_id to star (gold ★ + name), for the trait we came from.
# ---------------------------------------------------------------------------
phewas_plot <- function(df, levels, xlab = "Category", ylab = "-log10(P_mBATcombo)",
                        source = NULL, highlight = NULL, bonf = NULL, ymax = NULL) {
  ok <- !is.na(df$cat) & !is.na(df$P_mBATcombo)
  df <- df[ok, , drop = FALSE]
  if (!nrow(df)) {
    g <- plotly_empty(type = "scatter", mode = "markers")
    return(layout(g, annotations = list(list(
      text = "No associations", xref = "paper", yref = "paper",
      x = 0.5, y = 0.5, showarrow = FALSE,
      font = list(color = "#6c757d", size = 14)))))
  }
  
  levels <- levels[levels %in% df$cat]
  df$cat <- factor(df$cat, levels = levels)
  df     <- df[order(df$cat), , drop = FALSE]        # stable: keeps ascending-p within band
  df$y   <- pmin(-log10(df$P_mBATcombo), PVAL_CAP)
  
  # Equal-width slots: category k centres on integer k; its points spread across
  # [k-0.45, k+0.45] in ascending-p (rank) order.
  counts <- as.integer(table(df$cat))               # present levels, in order
  slot   <- rep(seq_along(counts), counts)          # 1..K per row
  idx    <- sequence(counts)                         # 1..n_i within each band
  n_i    <- rep(counts, counts)
  frac   <- ifelse(n_i == 1, 0.5, (idx - 1) / (n_i - 1))
  df$x   <- slot + (frac - 0.5) * 0.9
  
  # One rainbow hue per category (swept across the spectrum in level order).
  pal      <- grDevices::rainbow(length(counts), end = 0.85, v = 0.9)
  df$color <- pal[as.integer(df$cat)]
  
  # Evenly-spaced ticks at each slot centre.
  tlab <- df$label[match(levels(df$cat), as.character(df$cat))]
  
  df$hover <- paste0(
    "<b>", df$trait_name, "</b><br>",
    "P_mBATcombo: ", formatC(df$P_mBATcombo, format = "e", digits = 3), "<br>",
    "P_mBAT: ",      formatC(df$P_mBAT,      format = "e", digits = 3), "<br>",
    "P_fastBAT: ",   formatC(df$P_fastBAT,   format = "e", digits = 3)
  )
  
  bonf_y <- if (!is.null(bonf)) -log10(bonf) else NULL
  if (is.null(ymax)) ymax <- max(c(df$y, bonf_y), na.rm = TRUE) * 1.05
  
  g <- plot_ly(data = df, x = ~x, y = ~y, source = source, customdata = ~trait_id,
               type = "scattergl", mode = "markers",
               marker = list(color = df$color, size = 6),
               text = ~hover, hoverinfo = "text")
  g <- layout(g,
              xaxis = list(title = xlab, tickvals = seq_along(counts),
                           ticktext = as.character(tlab), tickangle = -45,
                           tickfont = list(size = 9)),
              yaxis = list(title = ylab, range = c(0, ymax)),
              showlegend = FALSE,
              margin = list(b = 140),
              shapes = if (!is.null(bonf_y)) list(list(
                type = "line", xref = "paper", x0 = 0, x1 = 1,
                yref = "y", y0 = bonf_y, y1 = bonf_y,
                line = list(color = "red", dash = "dash", width = 1))),
              annotations = if (!is.null(bonf_y)) list(list(
                xref = "paper", x = 1, xanchor = "right",
                yref = "y", y = bonf_y, yanchor = "bottom",
                text = paste0("Bonferroni p = ", formatC(bonf, format = "e", digits = 2)),
                showarrow = FALSE, font = list(color = "red", size = 11))))
  
  # Star the trait we came from (if it's in this plot's subset).
  if (!is.null(highlight) && nzchar(highlight)) {
    hi <- which(df$trait_id == highlight)
    if (length(hi) == 1) {
      g <- add_trace(g,
                     x = df$x[hi], y = df$y[hi], customdata = df$trait_id[hi],
                     type = "scattergl", mode = "markers+text", inherit = FALSE,
                     marker = list(color = "gold", size = 14, symbol = "star",
                                   line = list(color = "black", width = 1.5)),
                     text = df$trait_name[hi], textposition = "top center",
                     textfont = list(color = "black", size = 12),
                     hovertext = df$hover[hi], hoverinfo = "text",
                     showlegend = FALSE)
    }
  }
  
  g <- event_register(g, "plotly_click")
  g
}