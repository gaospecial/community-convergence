#' ---
#' title: "functions.R"
#' author: "Chun-Hui Gao"
#' note: "This is a collection of useful functions that I usually used"
#' ---
#'

# 在主成分分析图中添加椭圆形虚线框表示分组
# modified from FigureYa38 authored by Xu, Zhou-Geng
add_ellipase <- function(p, x="PC1", y="PC2", group="group",
                         ellipase_pro = 0.95,
                         linetype="dashed",
                         colour = "black",
                         lwd = 1,...){
  obs <- p$data[,c(x, y, group)]
  colnames(obs) <- c("x", "y", "group")
  ellipse_pro <- ellipase_pro
  theta <- c(seq(-pi, pi, length = 50), seq(pi, -pi, length = 50))
  circle <- cbind(cos(theta), sin(theta))
  ell <- plyr::ddply(obs, 'group', function(x) {
    if(nrow(x) <= 2) {
      return(NULL)
    }
    sigma <- var(cbind(x$x, x$y))
    mu <- c(mean(x$x), mean(x$y))
    ed <- sqrt(qchisq(ellipse_pro, df = 2))
    data.frame(sweep(circle %*% chol(sigma) * ed, 2, mu, FUN = '+'))
    })
  names(ell)[2:3] <- c('x', 'y')
  
  ell <- plyr::ddply(ell, plyr::.(group) , function(x) x[chull(x$x, x$y), ])
  p <- p + geom_polygon(data = ell, 
                        aes(x=x,y=y,group = group), 
                        colour = colour,
                        linetype=linetype,
                        lwd =lwd,
                        inherit.aes = F,
                        ...)
  return(p)
}

#' 生成指定物种的 count matrix 和 coldata
#'
#' @param ht_counts  
#' @param org: c("EC","PP")  
#'
#' @return a list of count matrix and coldata for DESeq2  
#' @export  
#'
#' @examples  
myDESeqMatrix <- function(ht_counts=ht_counts, org=NULL) {
  require(dplyr)
  require(tidyr)
  require(tibble)
  # 过滤和创建 matrix
  if (!is.null(org)){
    ht_counts <- ht_counts %>% filter(organism==org) %>%
      filter(ratio0 != ifelse(org=="EC","none","all"))
  }
  
  count_data <- ht_counts %>% dplyr::select(gene,count,sample_id) %>%
    group_by(sample_id) %>%
    spread(sample_id,count) %>% 
    as.data.frame() %>%
    column_to_rownames(var = "gene")
  col_data <- ht_counts %>% dplyr::select(sample_id,ratio0,time,group) %>% 
    arrange(sample_id) %>% distinct()
  return(list("count_data"=count_data,"col_data"=col_data))
}


#' 更好的科学计数法
#'
#' @param l
#'
#' @return
#' @export
#'
#' @examples
fancy_scientific <- function(l) {
  # human-friendly scitific notations
  # 
  # Turn in to character string in scientific notation
  # 
  # Usage:
  #   scale_y_continuous(labels=fancy_scientific,limits = c(2e-4,2e-3))
  
  l <- format(l, scientific = TRUE)
  # quote the part before the exponent to keep all the digits
  l <- gsub("^(.*)e", "'\\1'e", l)
  # turn the 'e+' into plotmath format
  l <- gsub("e", "%*%10^", l)
  # remove +
  l <- gsub("\\+","",l)
  # return this as an expression
  parse(text=l)
}



#' PCA 作图
#'
#' @param vsd: a matrix
#' @param show.label: if TRUE, sample name will be showed
#' @param return_data: if TRUE, return a dataframe of PCA results,
#'                     else, return the PCA plot.
#'
#' @return
#' @export
#'
#' @examples
myPlotPCA <- function(object,
                      intgroup = c("time","ratio0"), 
                      show.label = FALSE,
                      return_data = FALSE) {
  require(dplyr,quietly = T)
  require(DESeq2,quietly = T)
  require(vegan,quietly = T)
  require(ggrepel)
  # 进行 PCA 分析并返回 figure
  pca <- rda(t(assay(object))) 
  
  # 计算解释度
  percent_var <- pca$CA$eig/pca$tot.chi  # rda() 的结果中信息比较完整
  
  if (!all(intgroup %in% names(colData(object)))) {
    stop("the argument 'intgroup' should specify columns of colData(dds)")
  }
  intgroup_df <- as.data.frame(colData(object)[, intgroup, 
                                               drop = FALSE]) %>%
    tibble::rownames_to_column(var = "sample_id")
  
  # 提取数据和作图
  df <- scores(pca)$sites %>% 
    as.data.frame() %>%
    tibble::rownames_to_column(var="sample_id") %>%
    left_join(intgroup_df,by="sample_id") %>% 
    mutate(time=factor(time,levels = sort(unique(as.numeric(time)))))
  
  if (return_data){
    attr(df, "percentvar") <- percent_var
    return(df)
  } 
  
  mapping <- aes(PC1, PC2, color=ratio0, label=sample_id)
  
  p <- ggplot(df,mapping) +
    geom_point(size=2)  +
    xlab(paste0("PC1: ", round(percent_var[1] * 100), "% variance")) + 
    ylab(paste0("PC2: ", round(percent_var[2] * 100), "% variance"))
  
  if (show.label) {
    return(p + geom_text_repel(show.legend = F))
  } else {
    return(p)
  }
  
}

color_by_time <- function(name="Time(h)"){
  ggsci::scale_color_npg(name=name)
}

shape_by_ratio0 <- function(name="Initial Ratio\n(EC/PP)",...){
  scale_shape_manual(values = c(none=1, less=6, equal=5,more=2,all=0),
                     name = name,...)
} 

#' 根据样本距离矩阵绘制聚类树
#'
#' @param vsd: transformed data matrix
#' @param method: method used/passed to dist() function
#'                the distance measure to be used. 
#'                This must be one of "euclidean", "maximum", 
#'                "manhattan", "canberra", "binary" or "minkowski". 
#'                Any unambiguous substring can be given.
#' @param mapping: 
#'
#' @return ggtree object
#' @export
#'
#' @examples
myPlotTree <- function(object, method = "euclidean",
                       mapping=aes(shape = ratio0, color = factor(time)),...
                       ){
  require(ggtree,quietly = T)
  require(ape,quietly = T)
  # 计算样本间距
  sample_dists <- dist(t(assay(object)),method = method)
  
  # 绘制聚类图
  tree <- hclust(sample_dists) %>%   
    as.phylo()  # ape 可将 hclust 结果转变为 tree
  
  
  
  # 树 + 注释
  tr <- ggtree(tree) %<+% data.frame(colData(object))
  
  tr <- tr + geom_tiplab(mapping,offset = 3,show.legend = F) +
    geom_tippoint(mapping, size = 3) +
    # 这是让图例中的形状变大的办法
    # ref: https://stackoverflow.com/questions/15059093/ggplot2-adjust-the-symbol-size-in-legends
    guides(shape = guide_legend(override.aes = list(size = 3)),
           color = guide_legend(override.aes = list(size = 3))) +
    theme(legend.position = "right")
  xmax <- max(tr$data$x) * 1.1
  tr + xlim(NA,xmax)
}

#' 依据 comparison 的分组信息从 dds 中获取结果,并进行富集分析,返回dotplot
#'
#' @param dds 
#' @param comparison 
#' @param lfcThreshold 
#' @param p.adjusted 
#'
#' @return A list of DEG
#' @export
#'
#' @examples
myDEG_Results <- function(dds = NULL, comparison = NULL, 
                          lfcThreshold = 1,p.adjusted=0.05,
                          filtered = TRUE) {
  require(tibble,quietly = T)
  require(dplyr,quietly = T)
  results <- lapply(comparison, function(x){
    results(dds, contrast = c("group",x),
            lfcThreshold = lfcThreshold,
            alpha = p.adjusted)
  })
  names(results) <- sapply(comparison,function(x)paste(x,collapse = "_vs_"))
  for (i in seq_along(results)){
    results[[i]] %<>%
      as.data.frame() %>%
      rownames_to_column(var="gene") %>%
      as_tibble() %>%
      mutate(comparison = names(results[i])) 
  }
  if (!filtered) return(results)
  for (i in seq_along(results)){
    results[[i]] %<>%
      filter(padj < p.adjusted) %>%
      dplyr::select(gene,log2FoldChange,padj,comparison) %>%
      mutate(expression = ifelse(log2FoldChange>0,"up","dn")) 
  }
  return(results)
}


#' Get all named results of dds
#'
#' @param dds a DESeq object
#' @param ... passed to results()
#'
#' @return A dataframe with name column 
#' @export
#'
#' @examples
getAllNamedResults <- function(dds,...) {
  results_names <- resultsNames(dds)[-1]
  all.results <- vector("list",length(results_names))
  for (i in seq_along(results_names)){
    res <- results(dds,name = results_names[[i]],...)
    res <- as.data.frame(res) %>% 
      tibble::rownames_to_column(var="gene") %>% 
      mutate(name=results_names[[i]])
    all.results[[i]] <- res
  }
  do.call("rbind",all.results)
}

#' KEGG enrichment dotplot for clustered enrichment results
#'
#' @param DEG list of DEG
#' @param group 
#' @param othergroup 
#' @param organism 
#' @param ... passed to dotplot()
#'
#' @return
#' @export
#'
#' @examples
myKEGG_dotplot <- function(DEG, 
                           group="time", othergroup="ratio0", 
                           organism="eco",
                           showCategory=NULL,
                           ...) {
  require(clusterProfiler,quietly = T)
  require(ggplot2,quietly = T)
  DEG_df <- do.call("rbind",DEG) %>% 
    tidyr::separate(comparison,c("ratio0","time"),sep="_",extra="drop") %>% 
    mutate(time=as.numeric(gsub("h","",time))) %>%
    as_tibble()
  
  formula_res <- compareCluster(as.formula(paste("gene","~", group, "+", othergroup)),
                                data=DEG_df,fun="enrichMKEGG",organism=organism)

  p <- dotplot(formula_res,x=as.formula(paste("~",group)),showCategory=showCategory,...) + 
    facet_grid(as.formula(paste("~",othergroup)))
  
  #class(p$data[[group]]) <- class(DEG_df[[group]])
  #class(p$data[[othergroup]]) <- class(DEG_df[[othergroup]])
  
  p$data$time <- factor(p$data$time,levels = c(0,4,8,24))
  p$data$ratio0 <- factor(p$data$ratio0,levels=c("less","equal","more"))
  
  return(p)
}

myGO_dotplot <- function(DEG, group="time", othergroup="ratio0", 
                         OrgDb=org.EcK12.eg.db,
                         ont = "CC",
                         organism = "eco",
                         filterGO = TRUE,
                         simplifyGO = TRUE,
                         GOlevel = 4,
                         showCategory=NULL,
                         ...) {
  require(clusterProfiler,quietly = T)
  require(ggplot2,quietly = T)
  DEG_df <- do.call("rbind",DEG) %>% 
    tidyr::separate(comparison,c("ratio0","time"),sep="_",extra="drop") %>% 
    mutate(time=as.numeric(gsub("h","",time))) %>%
    as_tibble()
  kegg2ncbi <- bitr_kegg(DEG_df$gene,"kegg","ncbi-geneid",organism = organism)
  colnames(kegg2ncbi) <- c("gene","ncbi_id")
  DEG_df <- left_join(DEG_df,kegg2ncbi,by="gene")
  formula_res <- compareCluster(as.formula(paste("ncbi_id","~", group, "+", othergroup)),
                                data=DEG_df,
                                fun="enrichGO",
                                OrgDb=OrgDb,
                                ont=ont)
  
  # simplify(), gofilter() and dropGO() are useful in focus
  if (filterGO) formula_res <- gofilter(formula_res,level=GOlevel)
  if (simplifyGO) formula_res <- simplify(formula_res,...)
  
  p <- dotplot(formula_res,x=as.formula(paste("~",group)),showCategory=showCategory,...) + 
    facet_grid(as.formula(paste("~",othergroup)))
  
  #class(p$data[[group]]) <- class(DEG_df[[group]])
  #class(p$data[[othergroup]]) <- class(DEG_df[[othergroup]])
  
  p$data$time <- factor(p$data$time,levels = c(0,4,8,24))
  p$data$ratio0 <- factor(p$data$ratio0,levels=c("less","equal","more"))
  
  return(p)
  #return(formula_res) # for test
}

#' dotplot for enrichment result
#'
#' @param object 
#' @param x 
#' @param color 
#' @param size 
#' @param showCategory 
#' @param split 
#' @param title 
#' @param orderBy 
#' @param decreasing 
#'
#' @return ggplot
#' @export
#'
#' @examples
myDotPlot <- function(object, x="GeneRatio",color="p.adjust",
                      size="Count",
                      showCategory=10,split=NULL,
                      title="",orderBy="GeneRatio",decreasing=TRUE){
  colorBy <- match.arg(color,c("pvalue","p.adjust","qvalue","time","ratio0","expression"))
  df <- fortify(object,showCategory=showCategory,split=split) %>% 
    mutate(ratio0=factor(ratio0,levels = c("none","less","equal","more","all"))) %>%
    mutate(expression=factor(expression,levels = c("up","dn"))) %>%
    mutate(time=as.numeric(time))
  idx <- order(df[[orderBy]], decreasing = decreasing)
  df$Description <- factor(df$Description, levels=rev(unique(df$Description[idx])))
  ggplot(df,aes_string(x=x,y="Description",size=size,color=colorBy)) +
    geom_point() +
    ylab(NULL) + ggtitle(title)
}

#' Truncate string vector of ggplot axis label
#'
#' @param label 
#' @param maxLen 
#' @param dot 
#'
#' @return The first few words whose total length do not exceed maxLen
#' @export
#'
#' @examples
short_label <- function(label, maxLen = 50, dot = TRUE){
  l <- stringr::str_wrap(label, maxLen)
  
  if (dot) return(sub("\n.*"," ...", l))
  
  return(sub("\n.*","",l))
}


#' PCA plot with matrix
#'
#' @param m coldata is sites, rowdata is species.
#'
#' @return ggplot object
#' @export
#'
#' @examples
plotPCA.matrix <- function(m) {
  library(vegan)
  pca <- rda(t(m))
  percent_var <- pca$CA$eig/pca$tot.chi 
  df <- scores(pca)$sites  %>% 
    as.data.frame() %>% 
    tibble::rownames_to_column(var="sample") %>% 
    separate(sample,c("ratio0","rep"),sep="-")
  df$ratio0 <- factor(df$ratio0, levels = c("none","less","equal","more","all"))
  ggplot(df, aes(PC1,PC2,shape=ratio0,label=ratio0,color=ratio0))+
    geom_point(size=3,show.legend = F) +
    ggsci::scale_color_npg() +
    shape_by_ratio0(name="Intial Ratio") +
    xlab(paste0("PC1: ", round(percent_var[1] * 100), "% variance")) +
    ylab(paste0("PC2: ", round(percent_var[2] * 100), "% variance"))  + 
    directlabels::geom_dl(method="smart.grid")
}