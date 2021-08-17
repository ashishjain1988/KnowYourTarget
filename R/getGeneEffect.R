## To Supress Note
# utils::globalVariables(c(".","tempdir"))

#' Interactive plot to show the effects of genes on cancer proliferation
#' @description This function creates the interactive plot to show the effects of
#' genes on cancer proliferation using the DepMap CRISPR-Cas9 data.
#' different figures and plots explaining the MAF dataset.
#' @author Ashish Jain
#' @param genes A list of genes (Hugo symbol)
#' @param depMapVersion DepMap data version to be used
#' @param groupBySubtypes Group the cancer cell lines from cancer subtypes, default=FALSE
#' @export
#' @return A list containing plots for individual genes. For each gene, we have three plots showing
#' the DepMap CRISPR-Cas9 data in the form of histogram, boxplot, and combined plot.
#'
#' @examples
#' g<-getGeneEffectsInDepMap(genes=c("CHECK1))
#' g$CHECK1$combined
#'
#' @import depmap
#' @import ExperimentHub
#' @import ggplot2
#' @importFrom plotly ggplotly
#' @importFrom plotly subplot
#' @import tidyr

getGeneEffectsInDepMap<-function(genes=NULL,depMapVersion="",groupBySubtypes=FALSE){

  #genes<-c("CHEK1","MMP9","DCTN1","KLC4","HSP90AB1","CDKN1B","SOX9","CLN6","SYT8","NOS2")
  if (is.null(genes)) {
    stop("Need to enter the gene or gene list")
  }
  #require(depmap)
  #require(ExperimentHub)
  #require(tidyverse)
  #require(gridExtra)
  eh <- ExperimentHub::ExperimentHub()
  dep<-AnnotationHub::query(eh, "depmap")
  ##Get Crispr dependency Scores
  ##Version 21Q1
  # pedDepMapDataInfo<-read.table("~/myPART/DepMapPedSampleInfo.txt",sep = "\t",header = T)
  # pedDepMapDataInfo <- pedDepMapDataInfo[pedDepMapDataInfo$Class_for_Manuscript == "Pediatric",]
  crisprData<-eh[["EH5358"]]#depmap::depmap_crispr()
  metaData<-eh[["EH5362"]]#depmap::depmap_metadata()
  geneExpDataDepMap<-eh[["EH5360"]]
  #cnvDataDepMap<-eh[["EH5359"]]
  #mutationCallsDepMap<-eh[["EH5361"]]
  #t<-data.frame(crisprData %>% dplyr::filter(depmap_id %in% pedDepMapDataInfo$DepMap_ID))
  t<-data.frame(crisprData)
  plotsLists<-list()
  for(gene in genes)
  {
    #gene<-"CHEK1"
    #d$gene<-unlist(lapply(d$gene,FUN=function(x){return(unlist(str_split(x," "))[1])}))
    d <- t[t$gene_name %in% gene,]
    d1 <- merge(d,metaData,by="depmap_id",all.x=TRUE,all.y=FALSE) %>% dplyr::select(depmap_id,gene_name,primary_disease,subtype_disease,dependency,cell_line.x,cell_line_name)
    d1$subtype_disease[is.na(d1$subtype_disease)]<-d1$primary_disease[is.na(d1$subtype_disease)]
    medianDep<-d1 %>% group_by(primary_disease) %>% summarise(Mean=mean(dependency),Median=median(dependency))
    d2 <- merge(d1,geneExpDataDepMap[geneExpDataDepMap$gene_name %in% gene,],by="depmap_id",all.x=TRUE,all.y=FALSE) %>% dplyr::select(depmap_id,gene_name.x,primary_disease,subtype_disease,dependency,rna_expression,cell_line.x,cell_line_name)
    #print(paste0(genes,"-",round(min(medianDep$Median),digits = 2)))
    # d3<- merge(d2,mutationCallsDepMap[mutationCallsDepMap$gene_name %in% gene,],by="depmap_id",all.x=TRUE,all.y=FALSE) %>% dplyr::select(depmap_id,gene_name.x,primary_disease,subtype_disease,dependency,rna_expression,cell_line.x,cell_line_name,is_deleterious,is_tcga_hotspot,is_cosmic_hotspot)
    # d3$color[d3$is_deleterious==TRUE]<-"blue"
    # d3$color[d3$is_tcga_hotspot == TRUE]<-"red"
    # d3$color[d3$is_cosmic_hotspot == TRUE]<-"red"
    # d3$color[is.na(d3$color)]<-"black"
    d2$size<-ifelse(d2$rna_expression >=log(5+1),"high","low")
    # colorName<-unique(d3$color)
    g<-ggplot2::ggplot(d2,ggplot2::aes(x=primary_disease,y=dependency))+
      ggplot2::theme_classic(base_size = 10) +
      ggplot2::theme(text=ggplot2::element_text(face = "bold"),axis.text = ggplot2::element_text(size = 10),
                     axis.title = ggplot2::element_text(size = 15,face = "bold"),legend.background = ggplot2::element_rect(colour = "black")) +
      ggplot2::geom_boxplot()+
      ggplot2::geom_point(ggplot2::aes(text=sprintf("CellLine: %s<br>Subtype: %s<br>Primary: %s<br>Log2(TPM+1): %s<br>Dependency: %s", cell_line_name, subtype_disease,primary_disease,rna_expression,dependency),
                                       size=size))+
      ggplot2::scale_size_manual(values=c(2.2,1)) +
      #ggplot2::scale_color_manual(values=colorName) +
      ggplot2::guides(color="none") +
      ggplot2::theme(legend.position='none') +
      #ggplot2::geom_point(ggplot2::aes())+
      #ggplot2::geom_boxplot(outlier.colour="black", outlier.shape=16,outlier.size=0.5, notch=FALSE)+
      ggplot2::labs(x="", y = "Normalized Dependency Score") +
      ggplot2::geom_hline(yintercept=-1, linetype="dashed", color = "red") +
      #geom_dotplot(binaxis='y', stackdir='center',dotsize = 0.5) +
      ggplot2::coord_flip()
    ggly <- plotly::ggplotly(g,tooltip="text")

    g1<-ggplot2::ggplot(d1,aes(x=dependency)) +
      ggplot2::theme_classic(base_size = 10) +
      ggplot2::theme(text=element_text(face = "bold"),axis.text = element_text(size = 12),axis.title = element_text(size = 15,face = "bold"),legend.background = element_rect(colour = "black")) +
      ggplot2::labs(x="", y = "") +
      ggplot2::geom_histogram() +
      ggplot2::geom_vline(xintercept=-1, linetype="dashed", color = "red")
    #p <- cowplot::plot_grid(g1, g, align = "v",ncol = 1,rel_heights=c(1,2))
    p<-plotly::subplot(plotly::ggplotly(g1),ggly,nrows = 2,heights = c(0.25,0.75),titleX=TRUE)
    #save_plot(paste0("/Users/jaina13/myPART/WGSData/tumor-only-somatic-mafs/DepMap/Demap-AllCellLines-",round(min(medianDep$Median),digits = 2),"-CRISPR-",gene,".pdf"), p,base_height = 10,base_width = 8)
    plotsLists[[gene]]<- list(ggHistogram=g1,ggBoxplot=ggly,combined=p)
  }
  return(plotsLists)
}
