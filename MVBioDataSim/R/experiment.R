#' @title Simulation of a biological experiment
#' @description Simulates a whole multi-view microarray experiment.
#' @param regnet a regulatory network to be simulated
#' @param INITIAL a matrix containing initial values for each node of the network for each sample to be simulated
#' @param SIGNAL  a list containing a list of functions for each sample to be simulated
#' @param times timescale in which simulation will happen
#' @param sigmaEta standard deviation parameter of the model of technical noise
#' @param sigmaEps standard deviation parameter of the model of technical noise
simulateExperiment <- function(regnet, INITIAL, SIGNAL, times=1:100, 
                               sigmaEta=0.5, sigmaEps=3000)
{
  #input checking
  if(class(regnet) != 'RegulatoryNetwork')
    stop('regnet must be an object of type RegulatoryNetwork!')
  
  if(nrow(INITIAL) != length(Genes(regnet)) + length(Mirnas(regnet)))
    stop('INITIAL must have as many rows as the number of nodes of regnet!')
  if(is.null(rownames(INITIAL)))
    stop('INITIAL rownames must correspond to gene and miRNA names!')
  
  if(length(times) < 30 || max(times) < 80)
    warning('Low temporal resolution or not enough simulation time may degrade the quality of generated datasets.')
  
  if(sigmaEta < 0)
    stop('sigmaEta must be non negative!')
  
  if(sigmaEps < 0)
    stop('sigmaEps must be non negative!')
  
  #end input check
  n <- nrow(INITIAL)
  samples <- ncol(INITIAL)
  
  if(requireNamespace('tcltk', quietly=T))
  {  
    pb <- tcltk::tkProgressBar(title="Simulating Experiment...", min=0, max=samples, 
                      initial=0)
  }
  else
  {
    cat('Simulating Experiment...')
    pb <- txtProgressBar(min=0,max=samples,initial=0,style=3)
  }
  
  dataset <-list()
  dataset[["raw"]] <- matrix(0, nrow=n, ncol=samples)
  dataset[["log"]] <- matrix(0, nrow=n, ncol=samples)
  for(i in 1:samples)
  {
    pb$label(paste("Sample ", i, "...", sep=""))
    pb$up(i)
    ex <- .simulateNetwork(regnet=regnet, initial=INITIAL[, i], 
                           signal=SIGNAL[[i]], times=times, sigmaEta=sigmaEta, 
                           sigmaEps=sigmaEps)
    dataset[["raw"]][, i] <- ex[["raw"]]
    dataset[["log"]][, i] <- ex[["log"]]
  }
  close(pb)
  return(dataset)
}

.simulateNetwork <- function(regnet, initial, signal, times=1:100, 
                             sigmaEta=0.5, sigmaEps=3000)
{
  idx <- regnet$interactions != 0
  count <- sum(idx)
  
  hN <- Matrix(0, nrow=regnet$genes+regnet$mirnas, ncol=regnet$genes+regnet$mirnas)
  tN <- Matrix(0, nrow=regnet$genes+regnet$mirnas, ncol=regnet$genes+regnet$mirnas)
  gN <- Matrix(0, nrow=regnet$mirnas, ncol=regnet$genes)
  
  pN <- regnet$parameters$production + rnorm(regnet$genes+regnet$mirnas, 0, 0.05)
  dN <- regnet$parameters$degradation + rnorm(regnet$genes+regnet$mirnas, 0, 0.05)
  
  hN[idx] <- regnet$parameters$h[idx] + rnorm(count, 0, 0.05)
  tN[idx] <- regnet$parameters$theta[idx] + rnorm(count, 0, 0.05)
  
  mIdx <- idx[Mirnas(regnet), Genes(regnet)]
  mCount <- sum(mIdx)
  gN[mIdx] <- regnet$parameters$gamma[mIdx] + rnorm(mCount, 0, 0.05)
  
  dimnames(hN) <- dimnames(tN) <- dimnames(regnet$interactions)
  names(dN) <- names(pN) <- rownames(regnet$interactions)
  dimnames(gN) <- dimnames(regnet$parameters$gamma)
  
  pb <- txtProgressBar(min=min(times), max=max(times), initial=min(times), style=3)
  sampleEnv <- c(signal, list(pb=pb, rules=regnet$rules,hN=hN, tN=tN, gN=gN, pN=pN, dN=dN))
  
  ex <- rk4(initial, times, .regulation, sampleEnv)
  
  data <- .addNoise(ex[max(times), 2:ncol(ex)], sigmaEta=sigmaEta, 
                    sigmaEps=sigmaEps)
  return(data)
#   size <- regnet$genes + regnet$mirnas
#   genes <- signal$names
#   defPars <- regnet$parameters
#   rules <- regnet$rules
#   for(g in genes)
#   {
#     rules[[g]][[2]][[3]] <- substitute(signal[[name]](timePoint), list(name=g))
#   }
#   env <- defPars
#   m <- rep(0, size)
# 
#   env$production <- env$production + rnorm(size, m, sd=0.05)
#   env$check <- T
#   env$rules <- rules
#   env$signal <- signal
# 
#   ex <- rk4(initial, times, .regulation, env)
#   data <- .addNoise(ex[max(times), 2:ncol(ex)], sigmaEta=sigmaEta, 
#                     sigmaEps=sigmaEps)
#   return(data)
}

#Noise model described in 
#Rocke, D. M., & Durbin, B. (2001). A model for measurement error for gene 
#expression arrays. Journal of Computational Biology : A Journal of 
#Computational Molecular Cell Biology, 8(6), 557–69. 
#doi:10.1089/106652701753307485
#
#Durbin, B. P., Hardin, J. S., Hawkins, D. M., & Rocke, D. M. (2002). A 
#variance-stabilizing transformation for gene-expression microarray data. 
#Bioinformatics, 18(Suppl 1), S105–S110. 
#doi:10.1093/bioinformatics/18.suppl_1.S105
#
#NOTE: the parameter alpha of the original model is not included in this 
#implementation since we assume data is background-corrected.
.addNoise <- function(data, sigmaEta=0.5, sigmaEps=3000)
{
    size <- length(data)
    S_Eta_sq <- 2^(sigmaEta^2)*(2^(sigmaEta^2)-1)
    const <- sigmaEps^2/S_Eta_sq
    
    noise.mul   <- 2^rnorm(size, sd=sigmaEta)
    noise.add   <- rnorm(size, sd=sigmaEps)
    
    dataset <- list()
    dataset[["raw"]] <- (data*65535) * noise.mul + noise.add
    dataset[["log"]] <- log(dataset[["raw"]] + sqrt((dataset[["raw"]]^2)+const))
    return(dataset)
}

.zero <- function(t){0}
.one  <- function(t){1}

.regulation <- function(t, y, params)
{
  params$pb$up(t)
  output <- sapply(params$rules, eval, c(params, list(y=y, t=t)))
  #if(params$check)
  #{
    bad         <- (!is.finite(output) | (y<0))
    output[bad] <- 0
  output <- (1-y)*output
  #}
  return(list(output))
}