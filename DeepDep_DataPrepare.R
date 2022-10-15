#DeepDep data preparation


MB.gxp.deepmap<-cbind(row.names(MB.gxp),MB.gxp)
Prep4DeepDEP(exp.data = MB.gxp.deepmap, mut.data = NULL,
             meth.data = NULL,cna.data  = NULL,
             dep.data = NULL,mode = "Prediction",filename.out = "data_out" )

TC.gxp.deepmap<-cbind(row.names(TC.gxp),as.data.frame(TC.gxp))
Prep4DeepDEP(exp.data = TC.gxp.deepmap, mut.data = NULL,
             meth.data = NULL,cna.data  = NULL,
             dep.data = NULL,mode = "Prediction",filename.out = "data_out_tc" )




Prep4DeepDEP <- function(exp.data = NULL, mut.data = NULL,
                         meth.data = NULL,cna.data  = NULL,
                         dep.data = NULL, mode = c("Training","Prediction"),
                         filename.out = "data_out"){
  
  check.cellNames <- NULL
  cat("Mode",mode,"\n\n")
  path <- system.file("extdata/",package = "Prep4DeepDEP")
  
  #Check file
  if(sum(is.null(exp.data),is.null(mut.data),is.null(meth.data),is.null(cna.data)) == 4){
    stop(c("All genomic profiles are missing. Please provide at least one of mut.data, exp.data, meth.data, and cna.data."))
  }
  
  if(is.null(dep.data) & tolower(mode) == "prediction") {
    cat("dep.data is not provided, running with the default 1298 DepOIs...","\n")
    load(paste0(path,"default_dep_genes_1298.RData"))
  }
  
  if(is.null(dep.data) & tolower(mode) == "training") {
    cat("dep.data is not provided. Please provide gene dependency scores for the training mode...", "\n")
  }
  
  if(ncol(dep.data) == 1 & tolower(mode) == "training"){
    stop(c("Only one column detected in dep.data. Please provide gene dependency symbols and scores for the training mode."),call. = FALSE)
  }
  
  #Check gene symbol
  load(paste0(path,"gene_fingerprints_CGP.RData"))
  list.genes <- .CheckGeneSymbol(dep.data = dep.data,filename.out= filename.out)
  n <- nrow(list.genes)
  
  #Expression
  if(!is.null(exp.data)){
    if (!is.character(exp.data[1,1])) {
      stop("exp.data format error, please check!", call. = FALSE)
    }
    
    cat(c("Exp started..."),"\n")
    colnames(exp.data)[1] <- "Gene"
    ncell <- ncol(exp.data[,-1])
    if(is.null(check.cellNames)){
      check.cellNames <- colnames(exp.data[,-1])
    }else if(length(check.cellNames) != length( colnames(exp.data[,-1])) |
             sum(check.cellNames %in% colnames(exp.data[,-1])) != length(check.cellNames)){
      stop(c("Cell line names are inconsistent!"),call. = FALSE)
    }
    
    cat("Precessing",paste(length(check.cellNames),"cell lines..."),"\n")
    
    inputData <- exp.data[!duplicated(exp.data$Gene),]
    
    load(paste0(path,"ccle_exp_for_missing_value_6016.RData"))
    
    outputData <- merge(exp.index,inputData, by = "Gene", sort = FALSE, all.x = TRUE )
    Gene <- outputData$Gene
    rownames(outputData) <- outputData$Gene
    value_NA <-  rowSums(outputData[,-c(1,2)])
    
    cat(sum(is.na(value_NA)),"genes with NA values in exp.data. Substitute by mean values of CCLE.","\n")
    if(round((sum(is.na(value_NA))/nrow(outputData)),digits = 2)> 0.2){
      warning("NA found in >20% genes, please check if input format is correct!")
    }
    
    for (i in 1:nrow(outputData)) {
      if(is.na(value_NA[i])){
        outputData[i,is.na(outputData[i,])] <- outputData$Mean[i]
      }
      
    }
    outputData <- round(as.matrix(outputData[,-1]), digits = 4)
    outputData.final.exp  <- cbind(Gene,as.data.frame( outputData[,-1],stringsAsFactors = FALSE))
    outputData.final.exp  <- outputData.final.exp[,c("Gene", check.cellNames)]
    
    if(tolower(mode) == "prediction"){
      write.table(outputData.final.exp, file = paste(filename.out,"exp_prediction.txt",sep = "_"),
                  sep = "\t",row.names = FALSE,col.names = TRUE, quote = FALSE)
    }
    
    if (tolower(mode) == "training") {
      k = 2
      rep_col <- do.call("cbind", replicate(n,outputData.final.exp[,2], simplify = FALSE))
      colnames(rep_col) <- paste0("C1G",seq(1,n,1))
      if(ncol(outputData.final.exp) >= 3){
        for (i in 3:ncol(outputData.final.exp)) {
          rep_col.1 <- do.call("cbind", replicate(n,outputData.final.exp[,i],
                                                  simplify = FALSE))
          colnames(rep_col.1) <- paste0("C",k,"G",seq(1,n,1))
          k = k+1
          rep_col <- cbind(rep_col,rep_col.1)
        }
        
        rep_col <- cbind(rownames(outputData.final.exp),as.data.frame(rep_col,stringsAsFactors = FALSE))
        colnames(rep_col)[1] <- "Gene"
      }
      
      write.table(rep_col, file = paste(filename.out,"exp_training.txt",sep = "_"),
                  sep = "\t",row.names = FALSE,col.names = TRUE, quote = FALSE)
      
    }
    cat("Exp completed!","\n\n")
  }
  
  
  #Mutation
  if(!is.null(mut.data)){
    if (!is.character(mut.data[1,1])) {
      stop("mut.data format error, please check!", call. = FALSE)
    }
    cat(c("Mut started..."),"\n")
    
    colnames(mut.data)[1] <- "Gene"
    ncell <- ncol(exp.data[,-1])
    if(is.null(check.cellNames)){
      check.cellNames <- colnames(mut.data[,-1])
    }else if(length(check.cellNames) != length( colnames(mut.data[,-1])) |
             sum(check.cellNames %in% colnames(mut.data[,-1])) != length(check.cellNames)){
      stop(c("Cell line names are inconsistent!"),call. = FALSE)
    }
    cat("Precessing",paste(length(check.cellNames),"cell lines..."),"\n")
    
    inputData  <- mut.data[!duplicated(mut.data$Gene),]
    
    load(paste0(path,"ccle_mut_for_missing_value_4539.RData"))
    
    outputData <- merge(mut.index,inputData, by = "Gene", sort = FALSE, all.x = TRUE )
    Gene <- outputData$Gene
    rownames(outputData) <- outputData$Gene
    value_NA <-  rowSums(outputData[,-c(1,2)])
    cat(sum(is.na(value_NA)),"genes with NA values in mut.data. Substitute by median values of CCLE.","\n")
    if(round((sum(is.na(value_NA))/nrow(outputData)),digits = 2)> 0.2){
      warning("NA found in >20% genes, please check if input format is correct!")
    }
    
    for (i in 1:nrow(outputData)) {
      if(is.na(value_NA[i])){
        outputData[i,is.na(outputData[i,])] <- outputData$Median[i]
      }
      
    }
    
    outputData.final.mut  <- cbind(Gene, as.data.frame(outputData[,-c(1,2)]))
    outputData.final.mut  <- outputData.final.mut[,c("Gene",check.cellNames)]
    if(tolower(mode) == "prediction"){
      write.table(outputData.final.mut, file = paste(filename.out,"mut_prediction.txt",sep = "_"),
                  sep = "\t",row.names = FALSE,col.names = TRUE, quote = FALSE)
    }
    if (tolower(mode) == "training") {
      k = 2
      rep_col <- do.call("cbind", replicate(n,outputData.final.mut[,2]
                                            , simplify = FALSE))
      colnames(rep_col) <- paste0("C1G",seq(1,n,1))
      
      if(ncol(outputData.final.mut) >= 3){
        for (i in 3:ncol(outputData.final.mut)) {
          rep_col.1 <- do.call("cbind", replicate(n,outputData.final.mut[,i],
                                                  simplify = FALSE))
          colnames(rep_col.1) <- paste0("C",k,"G",seq(1,n,1))
          k = k+1
          rep_col <- cbind(rep_col,rep_col.1)
        }
        
        rep_col <- cbind(rownames(outputData.final.mut),as.data.frame(rep_col,stringsAsFactors = FALSE))
        colnames(rep_col)[1] <- "Gene"
      }
      
      write.table(rep_col, file = paste(filename.out,"mut_training.txt",sep = "_"),
                  sep = "\t",row.names = FALSE,col.names = TRUE, quote = FALSE)
    }
    cat("Mut completed!","\n\n")
  }
  
  
  #Methylation
  if(!is.null(meth.data)){
    if (!is.character(meth.data[1,1])) {
      stop("meth.data format error, please check!", call. = FALSE)
    }
    cat(c("Meth started..."),"\n")
    
    colnames(meth.data)[1] <- "Probe"
    ncell <- ncol(exp.data[,-1])
    if(is.null(check.cellNames)){
      check.cellNames <- colnames(meth.data[,-1])
    }else if(length(check.cellNames) != length( colnames(meth.data[,-1])) |
             sum(check.cellNames %in% colnames(meth.data[,-1])) != length(check.cellNames)){
      stop(c("Cell line names are inconsistent!"),call. = FALSE)
    }
    cat("Precessing",paste(length(check.cellNames),"cell lines..."),"\n")
    
    inputData  <- meth.data[!duplicated(meth.data$Probe),]
    
    load(paste0(path,"ccle_meth_for_missing_value_6617.RData"))
    
    outputData <- merge(meth.index,inputData, by = "Probe", sort = FALSE, all.x = TRUE )
    Probe <- outputData$Probe
    rownames(outputData) <- outputData$Probe
    value_NA <-  rowSums(outputData[,-c(1,2)])
    
    cat(sum(is.na(value_NA)),"genes with NA values in meth.data. Substitute by 0.","\n")
    if(round((sum(is.na(value_NA))/nrow(outputData)),digits = 2)> 0.2){
      warning("NA found in >20% genes, please check if input format is correct!")
    }
    
    for (i in 1:nrow(outputData)) {
      if(is.na(value_NA[i])){
        if(sum(is.na(outputData[i,])) == sum(ncol(outputData)-2)){
          # No values on all cell lines
          #stop(c(paste("Missing value",rownames(outputData)[i], sep = ":"),"\n"),call. = FALSE)
          outputData[i,is.na(outputData[i,])] <- 0
        }else{
          #No values on part of cell lines
          outputData[i,is.na(outputData[i,])] <- 0
        }
      }
    }
    
    outputData <- round(as.matrix(outputData[,-1]), digits = 4)
    outputData.final.meth  <- cbind(Probe, as.data.frame(outputData[,-1]))
    outputData.final.meth  <- outputData.final.meth[,c("Probe",check.cellNames)]
    
    if(tolower(mode) == "prediction"){
      write.table(outputData.final.meth, file = paste(filename.out,"meth_prediction.txt",sep = "_"),
                  sep = "\t",row.names = FALSE,col.names = TRUE, quote = FALSE)
    }
    
    if (tolower(mode) == "training") {
      k = 2
      rep_col <- do.call("cbind", replicate(n,outputData.final.meth[,2]
                                            , simplify = FALSE))
      colnames(rep_col) <- paste0("C1G",seq(1,n,1))
      
      if(ncol(outputData.final.meth) >= 3){
        for (i in 3:ncol(outputData.final.meth)) {
          rep_col.1 <- do.call("cbind", replicate(n,outputData.final.meth[,i],
                                                  simplify = FALSE))
          colnames(rep_col.1) <- paste0("C",k,"G",seq(1,n,1))
          k = k+1
          rep_col <- cbind(rep_col,rep_col.1)
        }
        
        rep_col <- cbind(rownames(outputData.final.meth),as.data.frame(rep_col,stringsAsFactors = FALSE))
        colnames(rep_col)[1] <- "Probe"
      }
      write.table(rep_col, file = paste(filename.out,"meth_training.txt",sep = "_"),
                  sep = "\t",row.names = FALSE,col.names = TRUE, quote = FALSE)
    }
    cat("Meth completed!","\n\n")
  }
  
  
  #Fingerprints
  if(!is.character(dep.data[1,1])){
    stop("dep.data format error, please check!", call. = FALSE)
  }else{
    cat(c("Fingerprint started..."),"\n")
    #load(paste0(path,"gene_fingerprints_CGP.RData"))
    
    colnames(dep.data)[1] <- "Gene"
    idx <- which(fingerprint[1,] %in% c("GeneSet",list.genes$Gene))
    outputData <- fingerprint[,idx]
    
    
    if(tolower(mode) == "prediction"){
      write.table(outputData, file = paste(filename.out,"fingerprint_prediction.txt",sep = "_"),
                  sep = "\t",row.names = FALSE,col.names = FALSE, quote = FALSE)
    }
    
    if(tolower(mode) == "training"){
      outputData.train <- cbind(outputData[,1],
                                do.call("cbind", replicate(ncell,outputData[,-1], simplify = FALSE)))
      outputData.train[1,] <- c("GeneSet",paste0(paste0("C",rep(seq(1,ncell,1),each = n)),"G",seq(1,n,1)))
      
      write.table(outputData.train, file = paste(filename.out,"fingerprint_training.txt",sep = "_"),
                  sep = "\t",row.names = FALSE,col.names = FALSE, quote = FALSE)
    }
    
    cat("Fingerprint completed!","\n\n")
  }
  
  #Dependency scores (training only)
  if (!is.null(dep.data) & tolower(mode) == "training") {
    if(is.null(check.cellNames)){
      check.cellNames <- colnames(dep.data[,-1])
    }else if(length(check.cellNames) != length( colnames(dep.data[,-1])) |
             sum(check.cellNames %in% colnames(dep.data[,-1])) != length(check.cellNames)){
      stop(c("Cell line names are inconsistent!"),call. = FALSE)
    }
    
    cat("Gene dependency scores (training mode) start...","\n")
    crispr.input <- dep.data[which(dep.data$Gene %in% list.genes$Gene),
                             which(colnames(dep.data) %in% c("Gene", check.cellNames))]
    crispr.output <- t(crispr.input[,2])
    colnames(crispr.output) <- paste0("C1G",seq(1,n,1))
    k=2
    for (i in 3:ncol(crispr.input)) {
      table <- t(crispr.input[,i])
      colnames(table) <- paste0("C",k,"G",seq(1,n,1))
      k = k+1
      crispr.output <- cbind(crispr.output,table)
    }
    
    crispr.output <- cbind("score",crispr.output)
    colnames(crispr.output)[1] <- "Dep_Score"
    write.table(crispr.output,file = paste(filename.out,"DepScore_training.txt",sep = "_"),
                sep = "\t",row.names = FALSE,col.names = TRUE, quote = FALSE)
    
    cat("Gene dependency scores (training) completed!","\n\n")
  }
  
  #Copy Number alteration
  if(!is.null(cna.data)){
    if (sum(colnames(cna.data) %in% c("CCLE_name","Chromosome","Start","End","Num_Probe", "Segment_Mean"))!=5) {
      stop("cna.data format error, please check!", call. = FALSE)
    }
    cat(c("CNA started..."),"\n")
    ncell <- length(unique(cna.data$CCLE_name))
    
    if(is.null(check.cellNames)){
      outputData.cna <- .PrepCNA(cna.original = cna.data,filename.out ,exportTable = FALSE)
    }else{
      idx <- which(cna.data$CCLE_name %in% check.cellNames)
      if(length(check.cellNames) != length(unique(cna.data$CCLE_name[idx])) |
         sum(check.cellNames %in% unique(cna.data$CCLE_name[idx])) != length(check.cellNames)){
        stop(c("Cell line names are inconsistent!"),call. = FALSE)
      }
      outputData.cna <- .PrepCNA(cna.original = cna.data[idx,], filename.out ,exportTable = FALSE)
      outputData.cna <- outputData.cna[,c("CNA",check.cellNames)]
    }
    
    if(tolower(mode) == "prediction"){
      colnames(outputData.cna)[1] <- "Bin"
      write.table(outputData.cna, file = paste(filename.out,"cna_prediction.txt",sep = "_"),
                  sep = "\t",row.names = FALSE,col.names = TRUE, quote = FALSE)
    }
    
    if (tolower(mode) == "training") {
      k = 2
      rep_col <- do.call("cbind", replicate(n,outputData.cna[,2]
                                            , simplify = FALSE))
      colnames(rep_col) <- paste0("C1G",seq(1,n,1))
      
      if(ncol(outputData.cna) >=3 ){
        for (i in 3:ncol(outputData.cna)) {
          rep_col.1 <- do.call("cbind", replicate(n,outputData.cna[,i],
                                                  simplify = FALSE))
          colnames(rep_col.1) <- paste0("C",k,"G",seq(1,n,1))
          k = k+1
          rep_col <- cbind(rep_col,rep_col.1)
        }
        rep_col <- cbind(outputData.cna$CNA,as.data.frame(rep_col,stringsAsFactors = FALSE))
        colnames(rep_col)[1]<- "Bin"
      }
      
      write.table(rep_col, file = paste(filename.out,"cna_training.txt",sep = "_"),
                  sep = "\t",row.names = FALSE,col.names = TRUE, quote = FALSE)
    }
    cat("CNA completed!","\n\n")
  }
}
