## To Supress Note
# utils::globalVariables(c(".","tempdir"))

#' Function to get the effects of target genes
#' @description This function creates an HTML file with the plots
#' using the NCI60 and DepMap data explaining the effects of the
#' input genes and drugs targetting them.
#' @author Ashish Jain
#' @param gene A gene to check the effect
#' @param outputFilePath Output path of the dashboard
#' @param outputFileName File name of the dashboard
#' @param dashboardTitle Title of the dashboard
#' @export
#' @return No return value, the dashboard html file is created
#' in the given output folder
#'
#' @examples
#' getGeneEffectsDashboard(gene="CHEK1")
#'
#' @importFrom i2dash i2dashboard
#' @importFrom i2dash datadir
#' @importFrom i2dash add_page
#' @importFrom i2dash add_component
#' @importFrom i2dash assemble
#' @importFrom DT datatable
#' @import rmarkdown

getGeneEffectsDashboard<-function(gene=NULL,outputFilePath=getwd(),outputFileName="KnowYourTarget.html",dashboardTitle="Know Your Target"){

  if (is.null(gene)) {
    stop("Need to enter the gene or gene list")
  }
  #require(i2dash)
  #require(plotly)
  topDrug <- getTopDrugs(genes = gene,noOfCores = 4)
  #topDrug1 <- topDrug$CHEK1[order(topDrug$CHEK1$Correlation),]
  topDrug1 <- topDrug[[gene]] %>% dplyr::filter(PValue <= 0.01) %>% arrange(Correlation,FDR)
  drugTable<- DT::datatable(topDrug1[,c(2,3,4,7)],rownames = FALSE,options = list(pageLength = 8))
  # drugTable <- plot_ly(
  #   type = 'table',
  #   #columnwidth = c(100, 100),
  #   #columnorder = c(0, 1),
  #   header = list(
  #     values = colnames(topDrug1),
  #     align = rep("center",ncol(topDrug1)),
  #     line = list(width = 1, color = 'black'),
  #     fill = list(rep("grey",ncol(topDrug1))),
  #     font = list(family = "Arial", size = 14, color = "white")
  #   ),
  #   cells = list(
  #     values = topDrug1,
  #     align = rep("center",ncol(topDrug1)),
  #     line = list(color = "black", width = 1),
  #     font = list(family = "Arial", size = 12, color = c("black"))
  #   ))
  nci60<-plotGeneDrugInteractionInNCI60(gene=gene,drug=topDrug1$drug[1])#"tolylquinone")
  depmapCrispr<-getGeneEffectsInDepMap(genes=gene)
  i2dash::i2dashboard(
    title = dashboardTitle,
    author = "",
    interactive = TRUE,
    theme = "yeti")->dashboard
  i2dash::datadir(dashboard) <- tempdir()
  dashboard %<>% i2dash::add_page(
      page = gene,
      title = gene,
      layout="2x2_grid",
      menu = NULL)
  dashboard %<>%
  i2dash::add_component(plotly::ggplotly(nci60[[paste0(topDrug1[1,]$drug,":",topDrug1[1,]$DrugNSC)]]),
                        page = gene,
                        title = "NCI60 Gene Expression VS Drug Activity") %>%
    i2dash::add_component(depmapCrispr[[gene]]$combined,
                          page = gene,
                          title = "Dependenacy Map using DepMap CRISPR-Cas9 data") %>%
  i2dash::add_component(drugTable,
                        page = gene,
                        title = "NCI60 Gene Expression VS Drug Activity Table")
  dashboard %>% i2dash::assemble(file = paste0(tempdir(),"MyDashboard.Rmd"), pages = gene)
  rmarkdown::render(paste0(tempdir(),"MyDashboard.Rmd"),
                    output_format="all", output_file=outputFileName,
                    output_dir = outputFilePath,
                    intermediates_dir = tempdir())
}
