## To Supress Note
# utils::globalVariables(c(".","tempdir"))

#' Function to get the top drugs using NCI60 data.
#' @description This function checkd for the top drugs targeting the
#' supplied genes using the NCI60 datasets.
#' @author Ashish Jain
#' @param genes A list of genes (Hugo symbol)
#' @param noOfCores No of cores used to run the function
#' @export
#' @return A data frame containing the correlation between the
#' gene expression and drug activity across all the drugs in
#' NCI60 data
#'
#' @examples
#' #result<-getTopDrugs(CHECK1"")
#' @import rcellminer
#' @import rcellminerData
#' @import doParallel
#' @importFrom foreach foreach
#' @import ggplot2
#' @importFrom lme4 lmer
#' @import dplyr

getTopDrugs<-function(genes=NULL,noOfCores=NULL){
  #require(rcellminer)
  #require(rcellminerData)
  if (is.null(genes)) {
    stop("Need to enter a gene or gene list")
  }
  if (is.null(noOfCores)) {
    stop("Need to enter a gene or gene list")
  }

  data("molData")
  geneExpMat <- exprs(molData[["exp"]])
  genesInNCI60Data<-intersect(row.names(geneExpMat),genes)
  sampleData<-getSampleData(molData)
  data("drugData")
  drugActMat <- exprs(getAct(drugData))
  drugAnnotDf <- as(featureData(getAct(drugData)), "data.frame")
  drugNames<-unique(drugAnnotDf$NAME)
  drugNames <- drugNames[drugNames != ""]
  #require(foreach)
  #require(doParallel)
  cores=parallel::detectCores()
  #setup parallel backend to use many processors
  if(is.null(noOfCores) || noOfCores > ncores)
  {
    noOfCores <- 1
  }

  cl <- parallel::makeCluster(noOfCores) #not to overload your computer cores[1]-4
  doParallel::registerDoParallel(cl)

  #row.names(geneExpMat)
  #genes<-c("FDFT1","LDLR","UPF1","ASCC3")
  #print(genes)
  finalCorrelation<-foreach::foreach(gene=genesInNCI60Data,.packages=c('tidyr','doParallel','foreach')) %dopar% {
    print(gene)
    finalMatrix <- foreach::foreach(x=drugNames,.packages='tidyr') %dopar% {
      drugPValue<-list()
      g<-geneExpMat[gene,]
      d1<-drugActMat[unlist(drugAnnotDf %>% dplyr::filter(NAME %in% x) %>% dplyr::select(NSC)), ,drop=FALSE]
      if(!is.null(nrow(d1)))
      {
        for(i in c(1:nrow(d1)))
        {
          corTest<-cor.test(g,d1[i,])
          l<-lm(g~d1[i,])
          l2<-lme4::lmer(g ~ d + (1 + d|t), data = data.frame(d=d1[i,],g=g,t=sampleData$TissueType))
          result<-c(paste0(x,":",row.names(d1[i,,drop=FALSE])),x,(corTest$estimate),(corTest$p.value),sigma(l),sigma(l2))
          names(result)<-c("DrugNSC","drug","Correlation","PValue","lmSigma","lmMixSigma")
          drugPValue[[i]]<-result
        }
      }
      return(do.call("rbind",drugPValue))
    }
    f<-data.frame(do.call("rbind",finalMatrix),stringsAsFactors = F)
    f$Correlation<-as.numeric(f$Correlation)
    f$PValue<-as.numeric(f$PValue)
    f$FDR<-p.adjust(f$PValue,method = "fdr")
    f<-f[order(f$FDR),]
    return(f)
  }
  names(finalCorrelation)<-genesInNCI60Data#row.names(geneExpMat)
  return(finalCorrelation)
}

#' Function returns the tissues in the NCI60 dataset
#' @description This function returns the tissues in the NCI60 dataset
#' @author Ashish Jain
#' @export
#' @return The list of tissues in the NCI60 dataset
#'
#' @examples

getAllTissuesInNCI60<-function(){
  #require(rcellminer)
  data("molData")
  geneExpMat <- exprs(molData[["exp"]])
  sampleData<-getSampleData(molData)
  return(sort(unique(sampleData$TissueType)))
}

#' Plot to show correlation between drug activity and gene expression
#' @description This function creates a plot to show the correlation between
#' the drug activity and gene expression
#' @param gene Gene name
#' @param drug Drug name
#'
#' @author Ashish Jain
#' @export
#' @return The list object containing the plot showing the correlation between
#' drug activity and gene expression
#'
#' @examples

plotGeneDrugInteractionInNCI60<-function(gene="CHEK1",drug="tolylquinone"){
  #require(rcellminer)
  #require(rcellminerData)
  g1<-NULL
  # Get the types of feature data in a MolData object.
  data("molData")
  geneExpMat <- exprs(molData[["exp"]])
  sampleData<-getSampleData(molData)
  data("drugData")
  drugActMat <- exprs(getAct(drugData))
  drugAnnotDf <- as(featureData(getAct(drugData)), "data.frame")
  d<-(drugAnnotDf %>% dplyr::filter(tolower(NAME) %in% tolower(drug)))
  listPlots<-list()
  if(nrow(d)>0)
  {
    for(r in c(1:nrow(d)))
    {
      #print(drug)
      g<-geneExpMat[gene,]
      d1<-drugActMat[d[r,1],]
      corTest<-cor.test(g,d1)
      g1<-ggplot2::ggplot(data.frame(GeneExp=g,DrugAct=d1,Tissue=stringr::str_to_title(sampleData$TissueType),
                                     CellLine=stringr::str_to_title(sampleData$Name),Cancer=stringr::str_to_title(sampleData$OncoTree2)),
                          ggplot2::aes(x=DrugAct,y=GeneExp)) +
        ggplot2::geom_point(ggplot2::aes(color=Tissue,label=Cancer,label2=CellLine)) +
        ggplot2::geom_smooth(method='lm',formula = y~x) +
        ggplot2::ggtitle(paste0(gene," VS ",d[r,2],"\nCor=",corTest$estimate,"\nP-Val:",corTest$p.value)) +
        ggplot2::theme_classic(base_size = 10) +
        ggplot2::labs(x="Drug Activity", y = "Gene Expression") +
        ggplot2::guides(color=ggplot2::guide_legend(title="")) +
        ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5),text=ggplot2::element_text(face = "bold"),
                       axis.text = ggplot2::element_text(size = 10),axis.title = ggplot2::element_text(size = 15,face = "bold"),
                       legend.background = ggplot2::element_rect(colour = "black"))
      #ggsave(paste0(parent,"/CancerSpecificGenes/WithDrugNames/NCI60-Plots/",tissue,"-",gene," VS ",d[r,2],"-",corTest$p.value,".pdf"),plot = g1)
      listPlots[[paste0(drug,":",row.names(d[r,]))]]<-g1
      # data.frame(GeneExp=g,DrugAct=d1,Tissue=stringr::str_to_title(sampleData$TissueType)) %>% filter(!is.na(GeneExp)) %>%
      #   plot_ly(x = ~DrugAct, y = ~GeneExp, mode = "markers") %>%
      #   add_markers(y = ~GeneExp) %>%
      #   add_trace(x = ~DrugAct, y = ~GeneExp, mode = "lines")
    }
  }
  return(listPlots)
}
