# 加载所需的包
library(bslib)
library(shiny)
library(shinythemes)
library(shinyjs)
library(plotly)
library(Seurat)
library(patchwork)
library(ggplot2)
library(readxl)
library(DT)
library(future)
library(bslib)
library(CATALYST)
library(openCyto)
library(flowWorkspace)
library(ggcyto)
library(flowCore)
library(impute)

options(future.globals.maxSize = 5 * 1024^3)  # 设置全局变量大小限制

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

if (!requireNamespace("glmGamPoi", quietly = TRUE))
  BiocManager::install("glmGamPoi")

source("ui/home_ui.R", local = TRUE)
source("ui/tab1_ui.R", local = TRUE)
source("ui/tab2_ui.R", local = TRUE)
source("ui/tab3_ui.R", local = TRUE)
source("ui/helptab_ui.R", local = TRUE)
source("server/tab1_server.R", local = TRUE)
source("server/tab2_server.R", local = TRUE)
source("server/tab3_server.R", local = TRUE)

options(shiny.maxRequestSize = 3000*1024^2)

# 通用 modal dialog 函数
showCompletionModal <- function(title = "Process Complete", message = "The process has finished successfully!") {
  showModal(modalDialog(
    title = title,
    message,
    easyClose = TRUE,
    footer = modalButton("OK")
  ))
}

ui <- page_navbar(
  # title = "SCPline",
  id = "nav",
  selected = "SCPline",  # c("SCPline", "Flow Cytometry", "Single Cell Multiomics", "MASS Workflow", "Help"),
  position = "static-top",
  collapsible = TRUE, # 是否允许折叠
  lang = "en",
  theme = bs_theme(version = 5, bootswatch = "yeti", font_scale = 1.1), 
  # 可以替换成：united，cerulean，journal，flatly，darkly，cosmo，cyborg，lumen，paper，sandstone，simplex，slate，spacelab，superhero，yeti，readable
  
  # bg = bs_get_variables(bs_theme(version = 5, bootswatch = "lumen"), "primary"),
  # nav_panel definitions
  nav_panel("SCPline", home_ui, icon = icon("home")),
  nav_panel("Flow Cytometry", tab1_ui, icon = icon("chart-simple")),
  nav_panel("Single Cell Multiomics", tab2_ui, icon = icon("chart-bar")),
  nav_panel("MASS Workflow", tab3_ui, icon = icon("chart-area")),
  nav_panel("Help", helptab_ui, icon = icon("question-circle"))
)

server <- function(input, output, session) {
  tab1_server(input, output, session)
  tab2_server(input, output, session)
  tab3_server(input, output, session)
}

shinyApp(ui = ui, server = server)
