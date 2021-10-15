# KnowYourTarget
Create an HTML Dashboard to know more about effects and drugs targeting a gene using NCI60 and DepMap data

**Requirement**

You need R version above 4.0.0 to run this application. Other dependencies are:

* `rcellminerData (>= 2.13)`
* `rcellminer (>= 2.12)`
* `rmarkdown (>= 2.9)`
* `depmap (>= 1.4)`
* `ggplot2 (>= 3.3.5)`
* `tidyr (>=1.1.3)`
* `dplyr (>= 1.0.2)`
* `knitr (>= 1.33)`
* `i2dash (>= 0.2.3)`
* `foreach (>= 1.5.1)`
* `doParallel (>= 1.0)`
* `ExperimentHub (>= 1.38.2)`
* `plotly (>= 4.9.4)`
* `DT (>= 0.18)`
* `stringr( >= 1.4)`
* `AnnotationHub (>=2.22)`

**How to install the R package**

**From Github**
* Install Dependencies
* `install.packages(c("dplyr","knitr","ggplot2","tidyr","foreach","doParallel","plotly","DT","stringr","rmarkdown","i2dash"))`
* `install.packages("BiocManager")`
* `BiocManager::install(c("rcellminerData","rcellminer","depmap","ExperimentHub","AnnotationHub"))`
* Now install the `devtools` package
* `install.packages("devtools")`
* `library(devtools)`
* Run command `install_github("ashishjain1988/KnowYourTarget")`

**More about the package**
 library(KnowYourTarget)

##Create a dashboard for a gene using NCI60 and DepMap data
getGeneEffectsDashboard(gene="CHEK1",outputFilePath=".",outputFileName="CHEK1-KnowYourTarget.html")

