# The helper function for plotting the plot. 

plot <- function(path, p, highlight = NULL){
  data <- readRDS(path)
  data$y <- -log10(data[[p]])

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
  
  graph <- plot_ly(data = data, x = ~cum_pos, y = ~y, 
                   type = "scatter", mode = "markers", marker = list(color = ~color), 
                   text = ~hover, hoverinfo = "text"
  )
  graph <- layout(graph,
                  xaxis = list(title = "Chromosome", tickvals = chr_mid, ticktext = names(chr_mid)),
                  yaxis = list(title = paste0("-log10(", p, ")")),
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
        type = "scatter", mode = "markers+text", inherit = FALSE,
        marker = list(color = "gold", size = 14, symbol = "star",
                      line = list(color = "black", width = 1.5)),
        text = data$Symbol[hi], textposition = "top center",
        textfont = list(color = "black", size = 12),
        hovertext = data$hover[hi], hoverinfo = "text",
        showlegend = FALSE)
    }
  }

  rm(data)
  return (graph)
}