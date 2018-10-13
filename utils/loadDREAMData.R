loadDreamDataSets <- function(root.dir = "~/Documents/uva/master/data/DREAM7/TrainSelect/"){
  require(qdapRegex)
  #List the available files
  file.names <- list.files(root.dir)
  
  #Ignore the readme files
  file.names <- file.names[!grepl(pattern = "README", x = file.names)]
  
  #Read in the data tables
  dream.data <- lapply(file.names, FUN = function(file.name){print(file.name);read.table(paste0(root.dir, file.name), header = TRUE, stringsAsFactors = FALSE)})
  
  #Get names for each file and set as names for the data list object produced above
  dream.names <- unlist(lapply(file.names, FUN = function(file.name){unlist(rm_between(file.name, left = "DREAM7_DrugSensitivity1_", right = ".txt", extract = TRUE))}))
  names(dream.data) <- dream.names
  
  #Process each data set (setting rownames and getting rid of meta data columns/rowname columns)
  dream.data <- lapply(dream.names, FUN = function(data.name){formatDreamData(data.name = data.name, dream.data = dream.data)})
  names(dream.data) <- dream.names
  
  return(dream.data)
}

formatDreamData <- function(data.name, dream.data){
  
  if(grepl("response", tolower(data.name))){
    data.set <- dream.data[[data.name]]
    rownames(data.set) <- data.set$CellLine
    data.set <- data.set[,-1]
    return(data.set)
  }else if (grepl("geneexpression", tolower(data.name))){
    data.set <- dream.data[[data.name]]
    rownames(data.set) <- data.set$HGNC_ID
    data.set <- data.set[,-1]
    return(data.set)
  }else if(grepl("methylation", tolower(data.name))){
    #Can also potentiall discretize the data here by thresholding the values at 0.2
    #the data as is represents the ratio of methylated to unmethylated for each probe sequence for each gene
    #thresholding would give us a binary indicator value indicating methylation of a gene
    data.set <- dream.data[[data.name]]
    rownames(data.set) <- with(data.set, paste0(Illumina_ID, ":", HGNC_ID))
    data.set <- data.set[, -c(1:4)]
    return(data.set)
  }else if(grepl("rnaseq", tolower(data.name))){
    data.set <- dream.data[[data.name]]
    rownames(data.set) <- with(data.set, paste0(Ensembl_ID, ":", HGNC_ID))
    data.set <- data.set[,-c(1:2)]
    return(data.set)
  }else if(grepl("rppa", tolower(data.name))){
    data.set <- dream.data[[data.name]]
    data.set <- data.set[data.set$FullyValidated %in% "Yes",]
    rownames(data.set) <- data.set$Antibody_ID
    data.set <- data.set[, -c(1:2)]
    return(data.set)
  }else if(grepl("snp", tolower(data.name))){
    data.set <- dream.data[[data.name]]
    rownames(data.set) <- with(data.set, paste0(EntrezID, ":", HGNC_ID))
    data.set <- data.set[,-c(1:2)]
    return(data.set)
  }
}