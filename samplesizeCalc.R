
# Do power calculation ----------------------------------------------------

run.sscalc <- function(z_a2, z_b, pi_0, treatment_effect, k, nr.percluster) {

  # pi_0 = proportion of population experiencing primary & second infection during trial at baseline
  # pi_a = proportion of population experiencing primary & second infection during trial in intervention group
  # k is between cluster variation coefficient
  
  pi_a <- pi_0 * (1-treatment_effect)
  
  clusters_perarm <- 2 + 
                    (z_a2 + z_b)^2 * 
                    ((pi_0 * (1-pi_0)/nr.percluster) + (pi_a * (1-pi_a)/nr.percluster) + k^2*(pi_0^2 + pi_a^2)) / (pi_0 - pi_a)^2 
  
  clusters_perarm <- ceiling(clusters_perarm)
  # print(clusters_perarm)
  return(clusters_perarm)
}

cl.perarm.inf <- run.sscalc(z_a2=1.96, z_b=0.84, pi_0=0.091, treatment_effect = 0.3, k=0.15, nr.percluster = 100)
cl.perarm.cases <- run.sscalc(z_a2=1.96, z_b=0.84, pi_0=0.015, treatment_effect = 0.3, k=0.15, nr.percluster = 20:200)
plot(cl.perarm.cases)
