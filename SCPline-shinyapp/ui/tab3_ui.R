tab3_ui <- page_sidebar(
  useShinyjs(),
  sidebar = sidebar(
    h4("MASS Workflow"),
    actionButton("btn_mass_overview3","Overview", class = "btn-primary btn-block", style = "margin-bottom: 1px;"),
    actionButton("btn_mass_upload", "Upload Data", class = "btn-primary btn-block", style = "margin-bottom: 1px;"),
    actionButton("btn_mass_preprocess", "Preprocess Data", class = "btn-primary btn-block", style = "margin-bottom: 1px;"),
    actionButton("btn_mass_pca", "PCA", class = "btn-primary btn-block", style = "margin-bottom: 1px;"),
    actionButton("btn_mass_umap", "UMAP", class = "btn-primary btn-block", style = "margin-bottom: 1px;"),
    actionButton("btn_mass_download", "Download Results", class = "btn-primary btn-block", style = "margin-bottom: 1px;"),
  ),
  div(
    id="tab3_main_content",
    div(
      id = "tab3_upload_panel",
      fluidRow(
        column(12,
               card(
                 card_header(h4("Load example data")),
                 card_body(
                   p("You can load example data to test the app."),
                   actionButton("btn_mass_load_example_data", "Run Example Data", class = "btn-primary btn-hover", style = "width: 20%")
                 )
               )
        ),
        column(12,
               card(
                 card_header(h4("Load User data")),
                 card_body(
                   fileInput("mass_peptide_data_path", label = "Upload your peptide data"),
                   fileInput("mass_metadata_path", label = "Upload your metadata"),
                   actionButton("btn_mass_load_user_data", "Load your Data", class = "btn-primary btn-hover", style = "width: 20%")
                 )
               )
        )
      )
    ),
    div(
      id = "tab3_preprocess_panel",
      fluidRow(
        column(12,
               card(
                 card_header(h4("Preprocess Data")),
                 card_body(
                   selectInput("mass_preprocess_method", "Preprocess Method", choices = c("sum", "mean", "median","medianPolish"), selected = "median"),
                   p("The data preprocessing includes the following steps:"),
                   tags$ul(
                     tags$li("Replace zeros with NA values for compatibility with downstream analysis."),
                     tags$li("Filter out data with more than 99% missing values to ensure data quality."),
                     tags$li("Count unique features at peptide and protein levels."),
                     tags$li("Aggregate peptide data to the protein level using various aggregation functions (sums, means, medians, and median polish)."),
                     tags$li("Normalize data by centering columns and rows with median and mean values."),
                     tags$li("Impute missing values using k-nearest neighbors (kNN) method.")
                   ),
                   actionButton("btn_mass_preprocess_data", "Preprocess Data", class = "btn-primary btn-hover", style = "width: 20%")
                 )
               ),
               card(
                 card_header(h4("Batch Effect Correction")),
                 card_body(
                   p("The batch effect correction includes the following steps:"),
                   tags$ul(
                     tags$li("Select a batch effect correction field."),
                     tags$li("Correct batch effects using the ComBat method.")
                   ),
                   selectInput("mass_batch_effect_correction_field", "Batch Effect Correction Field", choices = NULL, selected = NULL),
                   actionButton("btn_mass_batch_effect_correction", "Correct Batch Effects", class = "btn-primary btn-hover", style = "width: 20%")
                 )
               )
        )
      )
    ),
    div(
      id = "tab3_pca_panel",
      fluidRow(
        column(12,
               card(
                 card_header(h4("PCA")),
                 card_body(
                   p("The data analysis includes the following steps:"),
                   tags$ul(
                     tags$li("Convert scp data to Seurat object."),
                     tags$li("Perform principal component analysis (PCA) to visualize the data.")
                   ),
                   actionButton("btn_mass_pca_button", "Run PCA", class = "btn-primary btn-hover", style = "width: 20%")
                 )
               ),
               card(
                 card_header(h4("PCA Plot")),
                 card_body(
                   p("PCA is a dimensionality reduction technique that can be used for visualizing high-dimensional data in a low-dimensional space."),
                   plotlyOutput("mass_pca_plot", width = "400px", height = "400px")
                 )
               )
        )
      )
    ),
    div(
      id = "tab3_umap_panel",
      fluidRow(
        column(12,
               card(
                 card_header(h4("UMAP")),
                 card_body(
                   p("UMAP is a dimensionality reduction technique that can be used for visualizing high-dimensional data in a low-dimensional space."),
                   tags$ul(
                     tags$li("Perform clustering analysis to cluster the data."),
                     tags$li("Perform UMAP analysis to visualize the data.")
                   ),
                   actionButton("btn_mass_umap_button", "Run UMAP", class = "btn-primary btn-hover", style = "width: 20%")
                 )
               ),
               card(
                 card_header(h4("UMAP Plot")),
                 card_body(
                   p("UMAP is a dimensionality reduction technique that can be used for visualizing high-dimensional data in a low-dimensional space."),
                   selectInput("mass_umap_group_by", "Group by", choices = NULL, selected = NULL),
                   plotOutput("umap_plot", width = "400px", height = "400px")
                 )
               )
        )
      )
    ),
    div(
      id = "tab3_download_panel",
      h4("Download Results"),
      p("Please download your results here."),
      downloadButton("btn_mass_download_results", "Download", class = "btn-primary btn-hover"),
      img(
        src = "data/overview_img/footer.png", 
        style = "position: fixed; bottom: 10px; right: 10px; width: 500px; height: 120px; z-index: 1000;"
      )
    ),
    div(
      id = "tab3_overview3_panel",
      fluidRow(
        column(10,
               h4("Overview"),
               p("This section provides an overview of the MASS data analysis workflow."),
               img(src = "data/overview_img/tab3.overview.png", width = "800px")
        ),
        column(2,
               img(src = "data/overview_img/logo.png", width = "200px", height = "100px")
        )
      )
    )
  )
)