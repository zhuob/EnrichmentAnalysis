

test_sigPathway<- function(index, statistics, alternative = "either", type= "auto",
                           ranks.only = F, nsim=9999){
  alternative <- match.arg(alternative, c("mixed", "either", 
                                          "down", "up", "less", "greater", "two.sided"))
  if (alternative == "two.sided") 
    alternative <- "either"
  if (alternative == "less") 
    alternative <- "down"
  if (alternative == "greater") 
    alternative <- "up"
  type <- match.arg(tolower(type), c("auto", "t", "f"))
  allsamesign <- all(statistics >= 0) || all(statistics <= 
                                               0)
  if (type == "auto") {
    if (allsamesign) 
      type <- "f"
    else type <- "t"
  }
  if (type == "f" & alternative != "mixed") 
    stop("Only alternative=\"mixed\" is possible with F-like statistics.")
  if (alternative == "mixed") 
    statistics <- abs(statistics)
  if (alternative == "down") {
    statistics <- -statistics
    alternative <- "up"
  }
  if (ranks.only) {
    pvalues <- rankSumTestWithCorrelation(index = index, 
                                          statistics = statistics, df = Inf)
    p.value <- switch(alternative, down = pvalues["less"], 
                      up = pvalues["greater"], either = 2 * min(pvalues), 
                      mixed = pvalues["greater"])
  }
  else {
    ssel <- statistics[index]
    ssel <- ssel[!is.na(ssel)]
    nsel <- length(ssel)
    if (nsel == 0) 
      return(1)
    stat <- statistics[!is.na(statistics)]
    msel1 <- mean(ssel)
    msel2 <- mean(statistics[-index])
    msel <- msel1  - msel2
    
    if (alternative == "either") 
      posstat <- abs
    else posstat <- function(x) x
    msel <- posstat(msel)
    ntail <- 1
    permu.stat <- c()
    for (i in 1:nsim) {
      samp_id <- sample(1:500, nsel)
      permu.stat1 <- mean(stat[samp_id]) 
      permu.stat2 <- mean(stat[-samp_id]) 
      permu.stat[i] <- permu.stat1 - permu.stat2
      #      print(permu.stat[i])
      ntail <- ntail + (abs(permu.stat[i])> abs(msel))
    }
    p.value <- ntail/(nsim + 1)
    
  }
  return(list(p = p.value, obs.stat = msel, stat= permu.stat))
}

