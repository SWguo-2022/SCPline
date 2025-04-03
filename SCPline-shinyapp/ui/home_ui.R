home_ui <- fluidPage(
  id = "home_content",
  titlePanel("SCPline: An Interactive Shiny Framework for the Study of Single-Cell Proteomics Data Preprocessing"),
  
  # 左右结构：左边正文，右边图片
  fluidRow(
    column(6,
           div(
             id = "tab3_overview_panel",
             h4("Welcome to the Shiny Application of Single-Cell Proteomics Data Preprocessing!"),
             
             # 使用有序列表和段落分隔内容
             p("The Shiny application is designed to help researchers efficiently pre-process and analyze single-cell proteomics data. It is divided into the following three main modules:"),
             tags$ol(
               tags$li(
                 strong("Antibody-based method:"), 
                 " Used to process antibody labeling data from single-cell proteomic data. This module includes data cleaning, normalization, and visualization tools to help you quickly analyze antibody-labeled data."
               ),
               tags$li(
                 strong("Multi-omics method:"), 
                 " Based on Seurat, includes data screening, standardization, cluster analysis, and other data processing."
               ),
               tags$li(
                 strong("Mass spectrum-based approach:"), 
                 " Designed for processing mass spectrum-based single-cell proteomic data. Includes peptide polymerization method selection, standardization, missing value processing, and finally PCA and cluster analysis."
               )
             ),
             
             # 添加段落间的间距和强调部分
             p("The design goal of this application is to provide an intuitive interface and flexible data processing pipeline to help biological researchers who are not familiar with R programming better and more efficiently perform biology-related work.", style = "margin-top: 15px;"),
             
             # 使用 blockquote 提升可读性
             tags$blockquote(
               "Please use the navigation bar on the left to select the analysis module suitable for you and follow the steps in each module for data processing."
             ),
             
             # 加强呼吁操作的部分
             p("Click the buttons on the top to get started.", 
               style = "font-weight: bold; color: #0056b3; margin-top: 15px;")
           )
    ),
    column(6,
           img(src = "data/overview_img/homepage2.png", width = "100%", height = "30%"),
           img(src = "data/overview_img/homepage1.png", width = "100%", height = "30%"),
           img(src = "data/overview_img/homepage3.png", width = "100%", height = "30%")
    )
  )
)
