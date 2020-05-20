
mapToSymbol <- function(ensembls){
    ensembls <- gsub("\\..*", "", ensembls)
    symbols <- AnnotationDbi::mapIds(
    org.Hs.eg.db::org.Hs.eg.db,
    keys = ensembls,
    column ="SYMBOL",
    keytype ="ENSEMBL",
    multiVals ="first")
  
    for(i in 1:length(ensembls)){
        if(is.na(symbols[i])) 
            symbols[i] <- ensembls[i]
    }
  
    return(symbols)
}


mapToEnsembl <- function(symbols) {
  ensembls <- AnnotationDbi::mapIds(
    org.Hs.eg.db::org.Hs.eg.db,
    keys=symbols,
    column="ENSEMBL",
    keytype="SYMBOL",
    multiVals="first")

  ensemblIndicies <- which(gsub("\\..*", "", ensemblList) %in% ensembls)
  return(ensemblList[ensemblIndicies])
}
