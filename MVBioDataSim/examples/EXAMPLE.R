library(MVBioDataSim)

#create a regulatory network with default parameters
regnet <- RegulatoryNetwork(progress='gui')

#plot the regulatory network
plotRegulatoryNetwork(regnet)

#To simulate a simple dataset with two samples we must define the initial values
#for each gene and miRNA of the network. In addition we must specify a signal
#function for each controlling gene of the network for each samples

#We define random initial values for the two conditions
#initial values must have the same rownames of the genes and mirnas
INITIAL <- cbind(runif(regnet$genes+regnet$mirnas, 0,1),runif(regnet$genes+regnet$mirnas, 0,1))
rownames(INITIAL) <- rownames(regnet$interactions)

#Networks generated with default parameters have 5 controlling genes.
#NOTE: 	ALL SIGNALS 	AND EXPRESSION VALUES ARE BETWEEN 0 AND 1
#				ONLY FINAL EXPRESSION VALUES WILL COVER A MORE REALISTIC RANGE
#For the first sample we set to one the first three genes and knock down the remaining two
condition1 <- list(	s1=function(t){1},
										s2=function(t){1},
										s3=function(t){1},
										s4=function(t){0},
										s5=function(t){0})

#The second condition will be the opposite, except for the third gene that will be an oscillating signal
condition2 <- list(	s1=function(t){0},
										s2=function(t){0},
										s3=function(t){(1+sin(2 * pi * t / 50))/2},
										s4=function(t){1},
										s5=function(t){1})
										
SIGNALS <- list(condition1, condition2)

#Now the simulation will take place
ex <- simulateExperiment(regnet, INITIAL, SIGNALS)

#ex contains two matrices, raw correspond to raw expression data whereas log
#contains expression values with normalized variance as explained in
#Durbin, B. P., Hardin, J. S., Hawkins, D. M., & Rocke, D. M. (2002). 
#A variance-stabilizing transformation for gene-expression microarray data. 
#Bioinformatics, 18(Suppl 1), S105?S110. doi:10.1093/bioinformatics/18.suppl_1.S105