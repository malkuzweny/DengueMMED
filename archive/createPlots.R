# setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# source("dengue_working.R")

library(RColorBrewer)

# I4-I1 by age at equilibrium ---------------------------------------------

plotDF <- dengue_g_df
par(mar = c(5.1, 5.1, 4.1, 2.1))

font.size <- 2
lwd <- 4

# par(mar = c(bottom, left, top, right)) , alue for mar is c(5.1, 4.1, 4.1, 2.1)  

{
  cols <- brewer.pal(4, "Set1")
  
  plot( x=1:(max_age), y = plotDF[100,,"I1" ], type = "l", col=cols[1], lwd=lwd,
        ylim = c(0, 0.08), xlab="Age", ylab="% of population infected", cex.axis=font.size, cex.lab=font.size)
  
  lines( x = (1:max_age), y = plotDF[100,,"I2" ], col=cols[2], lwd=lwd)
  lines( x = (1:max_age), y = plotDF[100,,"I3" ], col=cols[3], lwd=lwd)
  lines( x = (1:max_age), y = plotDF[100,,"I4" ], col=cols[4], lwd=lwd)
  abline(v=4)
  abline(v=16)
  rect(xleft = 4, xright = 16, ybottom = par("usr")[3], ytop = par("usr")[4],
       border = NA, col = adjustcolor("black", alpha = 0.3))
  legend(x = c(25, 60), y = c(0.05, 0.08), legend=c("First infection", "Second infection", "Third infection", "Fourth infection"),
         col=c(cols[1], cols[2], cols[3], cols[4]), lty=1, cex=font.size, bg="white", bty="n", lwd=lwd)
  
}
 