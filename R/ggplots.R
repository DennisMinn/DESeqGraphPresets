#' Make Principle Component Analysis
#' 
#' Modified plotPCA from DESeq2 package that allows you to denote attributes/interest groups
#' with color, labels, or shapes.
#' @param object A \code{RangedSummarizedExperiment} object created by running regularized log transformation \code{DESeq2::rlog()} or
#' variance stabilizing transformation \code{DESeq2::vsd()} on a \code{DESeqDataSet} object.
#' @param color Desired intgroup/attribute from \code{RangedSummarizedExperiment::colData(obj)} you want denoted with color
#' @param label Desired intgroup/attribute from \code{RangedSummarizedExperiment::colData(obj)} you want denoted with labels
#' @param shape Desired intgroup/attribute from \code{RangedSummarizedExperiment::colData(obj)} you want denoted with shapes
#' @return Returns PCA either in your plots window, alternatively can be directed into a file with \code{pdf(file = "file_name")} and \code{dev.off()}
#' @examples 
#' > colnames(SummarizedExperiment::colData(rld)) 
#' [1] "replicate"  "sample"     "allele"     "sizeFactor"
#' 
#' ## Example 1 Label Sample Names
#' > mk_pca(rld)
#' 
#' ## Example 2 Highlight intgroup "replicate" by color
#' > mk_pca(rld, color = "replicate")
#' 
#' ## Example 3 Label samples by intgroup "replicate"
#' > mk_pca(rld, label = "replicate")
#' 
#' ## Example 4  Highlight all interest group by color, label, and shape
#' > mk_pca(rld, color = "replicate", label = "sample", shape = "allele")
#' 
#' ## Example 5 Exporting graph into a pdf 
#' > pdf("PCA.pdf")
#' > mk_pca(rld, color = "replicate", label = "sample", shape = "allele")
#' > dev.off()
#' @export 
mk_pca <- function(object, color, label, shape, ntop = 500){

    # Calculate the variance for each gene
    rv <- matrixStats::rowVars(SummarizedExperiment::assay(object))

    # Select the ntop genes by variance
    select <- order(rv, decreasing=TRUE)[seq_len(min(ntop, length(rv)))]

    # Perform a PCA on the data in assay(x) for the selected genes
    pca <- stats::prcomp(t(SummarizedExperiment::assay(object)[select,]))
    
    # Creates data frame for PCA
    df <- data.frame(PC1=pca$x[,1], PC2=pca$x[,2], name=colnames(object))

    # The contribution to the total variance for each component
    percentVar <- round(100 * pca$sdev^2 / sum( pca$sdev^2 ), digit = 2)
    attr(df, "percentVar") <- percentVar[1:2]

    # Handles provided arguments
    intgroup <- vector()
    if(!missing(color)){intgroup <- append(intgroup, color)}
    if(!missing(label)){intgroup <- append(intgroup, label)}
    if(!missing(shape)){intgroup <- append(intgroup, shape)}

    # Appends arguments to the PCA data frame
    for(group in intgroup){
        group_df <- data.frame(SummarizedExperiment::colData(object)[ , group, drop = FALSE])
        df <- cbind(df, group = group_df)
    }


    ggplot2::ggplot(df, ggplot2::aes_string("PC1", "PC2")) +
    ggplot2::xlab(paste0("PC1: ", percentVar[1],"% variance")) + #x-axis label 
    ggplot2::ylab(paste0("PC2: ", percentVar[2],"% variance")) + #y-axis label
    ggplot2::scale_colour_discrete() +                           # colour theme
    # Plots points and specifies shape and color if the arguments are provided 
    ggplot2::geom_point(
        size = 1.75, 
        ggplot2::aes_string(
            color = if(!missing(color)) color, 
            shape = if(!missing(shape)) shape)) +

    # Labels desired interest group if the argment is provided 
    {if(!missing(label)){
        ggrepel::geom_label_repel(
            size = 1.75,
            show.legend = FALSE,
            ggplot2::aes_string(
                label = label, 
                color = if(!missing(color)) color)) 
    }} +

    # Assign sample names to points if no other arguments are specified 
    {if(length(intgroup) == 0){
        ggrepel::geom_label_repel(
            size = 1.75,
            show.legend = FALSE,
            ggplot2::aes_string(label = df$name)) 
    }} 
}
#' Make Mean-Average Plot
#' 
#' Modified plotMA from DESeq2 package. Circles and labels either an input list of genes or genes with the lowest \code{n} p-adj values.
#' @param format Specifies which genes you would like to circle in \code{res}. \code{format="list"} will require that you input your desired list of gene into
#' \code{geneList}. \code{format="padj"} will require that you provide \code{res} and optionally \code{ntop} to specify the number of genes circled. 
#' @param res Contrast data produced by \code{DESeq2::res()} or \code{DESeq2::lfcShrink()}
#' @param pCutoff P-adj threshhold. Any values less than or equal to \code{pCutOff} will be marked in red, otherwise grey on the graph
#' @param geneList Input, desired list of genes
#' @param ntop Limits the number of genes included in \code{format = "padj"}. Additionally genes above \code{pCutoff} will not be highlighted
#' @param title Input title for graph
#' @param subtitle Input subtitle for graph
#' @param unit Labels graph to denote the comparison used for \code{res}
#' @param ylim Defines the y-axis boundaries 
#' @examples 
#' ## Example 1: Using padj
#' 
#' # MA-plot with the lowest p-adj genes circled (default 20 genes)
#' > mk_MAplot(format = "padj", res = res pCutOff = 0.05)
#' 
#' # Specifying lowest 10 p-adj genes 
#' > mk_mat(format = "padj", res = res, pCutoff = 0.05, ntop = 10)
#' 
#' # Setting the boundaries of the y-axis to [-5,5]
#' # Please make sure the ylim[1] is the lower bound while ylim[2] is the upper bound
#' > mk_mat(format = "padj", res = res, pCutoff = 0.05, ylim = c(-5,5))
#' 
#' ## Example 2: Using list
#' 
#' # Creates a matrix by subsetting the 'rld' with desired list of genes
#' > selectLab <- c("TSPAN6", "TNMD", "DPM1", "SCYL3", "C1orf112", "FGR")
#' > mat <- mk_mat(format = "list", res = res, geneList = selectLab)
#' 
#' ## Example 3: Exporting graph to pdf
#' > pdf("MAPlot.pdf")
#' > mk_MAplot(format = "padj", res = res pCutOff = 0.05)
#' > dev.off()
mk_MAplot <- function(format, res, pCutoff, geneList, ntop = 20, title, subtitle, unit, ylim){
    format <- match.arg(format, c("padj", "list"))
    if(missing(unit)){warning("No comparison specified for log2fold change")}
    if(!missing(ylim) && (length(ylim) != 2 || ylim[1] > ylim[2])){
        stop("Make sure ylim is a column vector with the first value representing \nthe lower bound the second value representing the upper bound")
    }

    # Maps ensembl IDs to gene symbols
    res$symbol <- mapToSymbol(rownames(res))
    # Creates a data frame and create "significant" column based on padj value
    res_df <- dplyr::mutate(as.data.frame(res), significant = padj <= pCutoff)

    # Handles NA values
    res_df$significant[is.na(res_df$significant)] <- FALSE
    res_df$log2FoldChange[is.na(res_df$log2FoldChange)] <- 0

    # Creates data frame for lowest p-adj genes
    if(format == "padj"){
        labels <- head(res_df[order(res$padj),],ntop)
        labels <- dplyr::filter(as.data.frame(labels), significant)
    } 
    # Creates data frame for input gene list
    if(format == "list"){
        labels <- which(res$symbol %in% geneList)
        labels <- res_df[labels,]
    }

    ggplot2::ggplot(res_df, ggplot2::aes(x=log2(baseMean), y=log2FoldChange, colour=significant)) +
    # Plots the points 
    ggplot2::geom_point(size = .5, cex = .35) +
    # Sets color 
    # Grey : Significant = FALSE 
    # Red  : Significant = TRUE
    ggplot2::scale_color_manual(values = c("grey", "red"))  +

    # Circle desired genes
    ggplot2::geom_point(
        data = labels, 
        size = .65, 
        cex = .35, 
        shape = 1, 
        color = "black",
        ggplot2::aes(
            x = log2(baseMean), 
            y = log2FoldChange)) +
    # Label desired genes
    ggrepel::geom_label_repel(
        data = labels, 
        color = "black", 
        size = 2, 
        segment.size = 0.1,
        ggplot2::aes(label=symbol)) +

    # Define Y-Axis Bounds
    {if(!missing(ylim))
    ggplot2::scale_y_continuous(
        breaks = seq(ylim[1],ylim[2],1), 
        limits = c(ylim[1],ylim[2]))} +

    # Subtitle and Y-Axis Label
    ggplot2::labs(
        font="Times New Roman",
        y = paste("log2(foldchange) ", if(!missing(unit)) unit), 
        subtitle = if(!missing(subtitle)) subtitle) +
    # Title
    ggplot2::ggtitle(label = if(!missing(title)) title) +
    # Removes Legend
    ggplot2::theme(legend.position = "none") 
}

#' Make Volcano Plot
#' 
#' Makes a volcano plot which circles and labels either an inputted list of gene or genes with the lowest \code{n} p-adj values.
#' Only genes that are both less than or equal to pCutoff and greater than or equal to the lfcCutoff will be denoted in red, otherwise grey.
#' @param format Specifies which genes you would like to circle in \code{res}. \code{format="list"} will require that you input your desired list of gene into
#' \code{geneList}. \code{format="padj"} will require that you provide \code{res} and optionally \code{ntop} to specify the number of genes circled. 
#' @param res Contrast data produced by \code{DESeq2::res()} or \code{DESeq2::lfcShrink()}
#' @param pCutOff P-adj threshhold. 
#' @param lfcCutoff Log2(foldchange) threshhold
#' @param geneList Input, desired list of genes
#' @param ntop Limits the number of genes included in \code{format = "padj"}. 
#' @param title title of graph
#' @param subtitle subtitle of graph
#' @param unit String to denote the comparison used for \code{res}
#' @param xlim Defines the x-axis boundaries 
#' @param ylim Defines the y-axis boundaries 
#' @examples 
#' ## Example 1: Using "p-adj"
#' 
#' # Volcano Plot with lowest p-adj genes, p-adj cutoff of 0.05, and log2(foldchange) of 1 (default 20 genes)
#' > mk_volcanoPlot(format = "padj", res = res, pCutoff = 0.05, lfcCutoff = 1) 
#' 
#' # Specifying lowest 10 p-adj genes
#' > mk_volcanoPlot(format = "padj", res = res, pCutoff = 0.05, lfcCutoff = 1, ntop = 10) 
#' 
#' # Specifying bounds for x-axis [-5,5] and y-axis [-9,10]
#' > mk_volcanoPlot(format = "padj", res = res, pCutoff = 0.05, lfcCutoff = 1, xlim = c(-5,5), ylim = c(-9, 10)) 
#' 
#' ## Example 2: Using "list"
#' 
#' > selectLab <- c("TSPAN6", "TNMD", "DPM1", "SCYL3", "C1orf112", "FGR")
#' > mk_volcanoPlot(format = "list", res = res, pCutoff = 0.05, lfcCutoff = 1, geneList = selectLab) 
#' 
#' ## Example 3: Exporting graph to pdf
#' > pdf("volcanoPlot.pdf")
#' > mk_volcanoPlot(format = "padj", res = res, pCutoff = 0.05, lfcCutoff = 1) 
#' > dev.off()
#' @export 
mk_volcanoPlot <- function(format, res, pCutoff, lfcCutoff, ntop = 20, geneList, unit, title, subtitle, xlim, ylim){
    format <- match.arg(format, c("padj", "list"))
    if(missing(unit)){warning("No comparison specified for log2fold change")}
    if(!missing(xlim) && (length(xlim) != 2 || xlim[1] > xlim[2])){
        stop("Make sure xlim is a column vector with the first value representing \nthe lower bound the second value representing the upper bound")
    }
    if(!missing(ylim) && (length(ylim) != 2 || ylim[1] > ylim[2])){
        stop("Make sure ylim is a column vector with the first value representing \nthe lower bound the second value representing the upper bound")
    }

    # Maps ensembl IDs to gene symbols
    res$symbol <- mapToSymbol(rownames(res))
    # Creates a data frame and create "significant" column based on pCutofff and lfcCutoff
    res_df <- dplyr::mutate(
        data.frame(res), 
        significant = padj <= pCutoff & abs(log2FoldChange) >= lfcCutoff)

    # Creates data frame for lowest p-adj genes
    if(format == "padj"){
        labels <- head(res_df[order(res$padj),],ntop)
        labels <- dplyr::filter(data.frame(labels), significant)
    } 

    # Creates data frame for input gene list
    if(format == "list"){
        labels <- which(res$symbol %in% geneList)
        labels <- res_df[labels,]
    }

    # Creates labels for x-axis and y-axis
    y_lab <- "-log10(p-adj)"
    x_lab <- paste("log2(foldchange) ", if(!missing(unit)) unit)

    ggplot2::ggplot(res_df, ggplot2::aes(x=log2FoldChange, y=-log10(padj), colour=significant)) +
    # Plots the points
    ggplot2::geom_point(size = 1, cex = 1) +
    # Sets color 
    # Grey : Significant = FALSE 
    # Red  : Significant = TRUE
    ggplot2::scale_color_manual(values = c("grey", "red"))  +
    
    #Circles desired genes
    ggplot2::geom_point(
        data = labels, 
        ggplot2::aes(x=log2FoldChange, 
        y = -log10(padj)), 
        size = 1, 
        cex = 1, 
        shape=1, 
        color="black") +
    #Labels desired genes
    ggrepel::geom_label_repel(
        data = labels, 
        ggplot2::aes(label=symbol),  
        size = 1.5, 
        segment.size = 0.1,
        color = "black") +

    #Draws P-adj Cutoff
    ggplot2::geom_hline(ggplot2::aes(yintercept = -log10(pCutoff)), linetype = "dashed") +
    #Draws Log2FoldChange Cutoff
    ggplot2::geom_vline(ggplot2::aes(xintercept = -lfcCutoff), linetype = "dashed") +
    ggplot2::geom_vline(ggplot2::aes(xintercept =  lfcCutoff), linetype = "dashed") +

    # Define Y-Axis Bounds
    {if(!missing(ylim))
    ggplot2::scale_y_continuous(
        breaks = seq(ylim[1],ylim[2],1), 
        limits = c(ylim[1],ylim[2]))} +

    # Define X-Axis Bounds
    {if(!missing(xlim))
    ggplot2::scale_x_continuous(
        breaks = if(!missing(xlim)) seq(xlim[1],xlim[2],1), 
        limits = if(!missing(xlim)) c(xlim[1],xlim[2]))} +

    # Subtitle, X-Axis, and Y-Axis Labels 
    ggplot2::labs(
        y = y_lab, 
        x = x_lab,
        subtitle = if(!missing(subtitle)) subtitle, 
        font="Times New Roman") +

    #Title
    ggplot2::ggtitle(label = if(!missing(title)) title) +
    #Removes legend
    ggplot2::theme(legend.position = "none")
}