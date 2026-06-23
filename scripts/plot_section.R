# The helper function for plotting the plot. 

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