
#' initDb - Initialization of the package so the NABEC Cortex network can be used with
#' CoExpNets
#'
#' @param mandatory If this parameter is `TRUE` then the networks will be added no matter whether they were already there.
#'
#' @return No value
#' @export
#'
#' @examples
initDb = function(mandatory=F){
  category = "nabec"
  the.dir = system.file("", category, package = "CoExpNABEC")
  tissue = getNABECTissues()
  
    CoExpNets::addNet(which.one=category,
           tissue=tissue,
           netfile=getNABECNet(tissue),
           ctfile=paste0(getNABECNet(tissue),".celltype.csv"),
           gofile=paste0(getNABECNet(tissue),"_gprofiler.csv"),
           exprdatafile=paste0(the.dir,"/",tissue,"_resids.rds"),
           overwrite=mandatory)
}


#' getNABECNet - Accessing a network object directly
#'
#' @param tissue One of the labels that can be obtained by calling getGTExTissues() method
#' to refer to a specific network within the package
#'
#' @return The RDS full file name to the network
#' @export
#'
#' @examples
getNABECNet = function(tissue=getNABECTissues()){
  stopifnot(tissue == getNABECTissues())
  the.dir = system.file("", "nabec", package = "CoExpNABEC")
  return(paste0(the.dir,"/netFCTX.30P.8.it.30.rds.gene.rds"))
}

#' Title Getting all NABEC available tissues
#'
#' @return
#' @export
#'
#' @examples
getNABECTissues = function(){
  return("Cortex")
}

prepareResiduals = function(){
  trasposeDataFrame = function(file.in,first.c.is.name=T){
    if(typeof(file.in) == "character")
      data.in = readRDS(file.in)
    else{
      data.in = file.in
      rm(file.in)
    }
    
    if(first.c.is.name){
      data.t = as.data.frame(cbind(apply(data.in[,-1],MARGIN=1,function(x){ return(as.numeric(x))})))
      colnames(data.t) = data.in[,1]
      rownames(data.t) = colnames(data.in[-1])
    }else{
      data.t = as.data.frame(cbind(apply(data.in,MARGIN=1,function(x){ return(as.numeric(x))})))
      colnames(data.t) = rownames(data.in)
      rownames(data.t) = colnames(data.in)	
    }
    return(data.t)
  }
  getUsableTranscriptsFromDict = function(dictf){
    if(typeof(dictf) == "character")
      dictionary = readRDS(dictf)
    else
      dictionary = dictf
    
    return(unlist(lapply(dictionary,function(x){
      mods = unique(x)
      ts = NULL
      for(mod in mods){
        t.in.mod = names(x)[x == mod]
        
        #Which transcript to use???
        #At the moment alphabetical order
        ts = c(ts,sort(t.in.mod)[1])
      }
      return(ts)
    })))
  }
  
  dictf="~/Dropbox/KCL/workspace/CoExpNABEC/inst/nabec/dictionary.1.0.7.rds"
  exprf = "~/Dropbox/KCL/workspace/CoExpNABEC/inst/nabec/autosomal.normalized.peer.adjusted"
  samplef = "~/Dropbox/KCL/workspace/CoExpNABEC/inst/nabec/nabec.postQC.core.cohort.list"
  expr.raw = read.table(exprf,stringsAsFactors=F)
  expr.data = trasposeDataFrame(expr.raw,F)	
  transcripts = getUsableTranscriptsFromDict(dictf)
  cat("We'll start from a dictionary of",length(transcripts),"\n")
  samples = read.table(samplef,stringsAsFactors=F)$V1
  samples = gsub("-",".",samples)
  expr.data = expr.data[rownames(expr.data) %in% samples,colnames(expr.data) %in% transcripts]
  cat("We'll create a network from",nrow(expr.data),"and",ncol(expr.data),"genes\n")
  saveRDS(expr.data,"~/Dropbox/KCL/workspace/CoExpNABEC/inst/nabec/nabec_resids.rds")
  
}
