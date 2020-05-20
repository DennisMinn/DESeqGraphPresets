#' Make Matrix
#'
#' Creates a matrix that highlights the upregulation and downregulation of either an input list of genes or genes with the lowest \code{n} p-adj values. 
#' Inspire by the BioConductor's course material "Differential expression, manipulation, and visualization of RNA-seq reads"
#' by Michael Love, Simon Anders, Wolfgang Huber.
#' @param format Specifies how to subset the matrix. \code{"list"} will produce a matrix with entries corresponding to the genes symbols
#' provide in \code{geneList}. \code{"padj"} will subset \code{object} based on the metadata \code{'padj'} provided by \code{res}.  
#' @param object A \code{RangedSummarizedExperiment} object created by running regularized log transformation \code{DESeq2::rlog()} or
#' variance stabilizing transformation \code{DESeq2::vsd()} on a \code{DESeqDataSet} object.
#' @param geneList Desired list of gene symbols 
#' @param res \code{DESeqResults} object \code{DESeq2::results()} containing the necessary metadata to subset \code{object}
#' @param pCutOff P-adj threshhold. Because Heatmaps do not indicate p-adj values, any genes above \code{pCutOff} will not be included in the matrix
#' @param nTop Limits the number of genes that are included in the matrix in \code{format = "padj"}, by default this value is set to 100.
#' @return Returns a matrix that highlights the the upregulation and downregulation of gene expression across samples; the values
#' within each entry are subtracted from the entry's mean values. That way, upregulation are indicated by positive values while downregulation are indicated by negative values. 
#' @examples 
#' ## Example 1: Using padj
#' 
#' # Creates a matrix by subsetting the 'rld' the with p-adj threshhold of 0.05 (defaults to 100 rows)
#' > mat <- mk_mat(format = "padj", object = rld, res = res, pCutoff = 0.05)
#' 
#' # Creates a matrix by subsetting the 'rld' with 20 lowest p-adj genes, granted the values are under 0.05
#' > mat <- mk_mat(format = "padj", object = rld, res = res, pCutoff = 0.05, nTop = 20))
#' 
#' ## Example 2: Using list
#' 
#' # Creates a matrix by subsetting the 'rld' with desired list of genes
#' > selectLab <- c("TSPAN6", "TNMD", "DPM1", "SCYL3", "C1orf112", "FGR")
#' > mat <- mk_mat(format = "list", object = rld, geneList = selectLab)
#' @export

mk_mat <- function(format, object, geneList, res, pCutoff, nTop = 100){
    
    if(missing(object)){stop("Missing RangedSummarizedExperiment object")}
    format <- match.arg(format, c("padj", "list"))

    if(format == "padj"){
        #Checks if pCutoff is valid
        stopifnot(pCutoff >= 0)

        #Reorders the RangedSummarizedExperiment based on pCutoff and then adds all entries 
        #that meet the threshholds into a matrix
        mat <- SummarizedExperiment::assay(object)[0, ]
        res_padj <- res[order(res$padj), 5, drop = FALSE]
        for(i in 1:n){
            if(res_padj[i,] >= pCutoff)
                break;
            ensemblID <- rownames(res_padj)[i]
            mat <- rbind(mat, SummarizedExperiment::assay(object)[rownames(object) %in% ensemblID, ])
            rownames(mat)[i] <- ensemblID
        }
    }

    if(format == "list"){
        #initial argument testing
        if(missing(geneList)){stop("Missing gene input")}
        
        #After converting the symbols into ensembl IDs, we find
        #the location of each gene in the RangedSummarizedExperiment and 
        #subset the values.
        ensembls <- mapToEnsembl(geneList)
        ensemblIndicies <- which(rownames(object) %in% ensembls)
        mat <- SummarizedExperiment::assay(object)[ensemblIndicies, ]
    }

    #To better highlight the upregulation and downregulation of gene expression
    #across samples, the values within each entry are subtracted from its mean.
    #That way, upregulation are indicated by positive values while downregulation are
    #indicated by negative values.
    mat <- mat - rowMeans(mat)

    #Mapping all ensembl IDs into symbols
    rownames(mat) <- mapToSymbol(rownames(mat))
    return(mat)
}

#' Make Data Frame
#' 
#' Creates the legend for the heatmap
#' @param attributes column vector containing the names of desired metadata
#' @param object A \code{RangedSummarizedExperiment} object created by running regularized log transformation \code{DESeq2::rlog()} or
#' variance stabilizing transformation \code{DESeq2::vsd()} on a \code{DESeqDataSet} object.
#' @return Returns a data frame that will be used as \code{df} argument in \code{mk_heatmap}
#' @examples
#' ## Example 1
#' 
#' # Creates a data frame using all metadata given from your design table in your 'object'
#' > df <- mk_df(attributes = c("replicate", "sample", "allele"), object = rld)
#' 
#' ## Example 2
#' 
#' # Creates a data frame using only two attributes give form your design table in your 'object'
#' > df <- mk_df(attributes = c("replicate", "sample"), object = rld)
#' 
#' @export
mk_df <- function(attributes, object){
    return(data.frame(SummarizedExperiment::colData(object)[ , attributes]))
}

#' Make Heatmap
#' 
#' Creates a heatmap inspire by the BioConductor's course material "Differential expression, manipulation, and visualization of RNA-seq reads"
#' by Michael Love, Simon Anders, Wolfgang Huber using \code{pheatmap} package with preset values 
#' @param mat matrix containing the expressions of desired genes, look at \code{mk_mat} for more detail
#' @param df data frame formatted from your design table, look at \code{mk_df} for more detail
#' @param sampleNames column vector containing the names of your samples
#' @param title optional argument containing the title of your graph
#' @return Returns a heatmap either in your plots window, alternatively can be directed into a file with \code{pdf(file = "file_name")} and \code{dev.off()}
#' @examples
#' > colnames(SummarizedExperiment::colData(rld)) 
#' [1] "replicate"  "sample"     "allele"     "sizeFactor"
#' > sampleNames 
#' [1] "sample.1"  "sample.2"  "sample.3"  "sample.4"  "sample.5"  "sample.6"  "sample.7"  "sample.8"  "sample.9" 
#' [10] "sample.10" "sample.11" "sample.12" "sample.13" "sample.14" "sample.15" "sample.16" "sample.17" "sample.18"  
#' 
#' ## Example 1 Creating a heatmap
#' > mk_heatmap(mat = mat, df = df, sampleNames = sampleNames, title = "heatmap")
#' 
#' ## Example 2 Directs the output into a pdf
#' > pdf(file = "Heatmap.pdf")
#' > sampleNames <- paste("sample", 1:18, sep = ".")
#' > mk_heatmap(mat = mat, df = df, sampleNames, "heatmap")
#' > dev.off()
#' 
#' @export 
mk_heatmap <- function(mat, df, sampleNames, title){
    if(missing(title)){title <- ""}
    pheatmap::pheatmap(
    mat,                        #Matrix of count data 
    annotation_col=df,          #Creates legend
    main = title,               #Creates title
    labels_col = sampleNames,  #Sample Names
    cellheight = 4,             #Height of displayed data
    cellwidth = 8,              #Width of displayed data
    fontsize = 4)               #Font size of genes
}
