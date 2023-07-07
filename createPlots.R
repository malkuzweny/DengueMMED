#setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

#source("dengue_working.R")

library(RColorBrewer)

# I4-I1 by age at equilibrium ---------------------------------------------

plotDF <- dengue_g_df

  
{
  cols <- brewer.pal(4, "Set1")
  
  plot( x=1:(max_age), y = plotDF[100,,"I1" ], type = "l", col=cols[1],
        ylim = c(0, 0.08), xlab="Age", ylab="% of population infected")
  
  lines( x = (1:max_age), y = plotDF[100,,"I2" ], col=cols[2])
  lines( x = (1:max_age), y = plotDF[100,,"I3" ], col=cols[3])
  lines( x = (1:max_age), y = plotDF[100,,"I4" ], col=cols[4])
  abline(v=4)
  abline(v=16)
  rect(xleft = 4, xright = 16, ybottom = par("usr")[3], ytop = par("usr")[4],
       border = NA, col = adjustcolor("black", alpha = 0.3))
  legend("topright", legend=c("First infection", "Second infection", "Third infection", "Fourth infection"),
         col=c(cols[1], cols[2], cols[3], cols[4]), lty=1, cex=0.8, bg="white")
  
}
