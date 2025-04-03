tab2_ui <- page_sidebar(
  useShinyjs(),
  sidebar = sidebar(
    h4("Multi-omics"),
    actionButton("btn_tab2_overview", "Overview", class = "btn btn-primary btn-block"),
    actionButton("btn_tab2_upload_data", "Upload your data", class = "btn btn-primary btn-block"),
    actionButton("btn_tab2_norm_feature", "Normalization & Feature Selection", class = "btn btn-primary btn-block", disabled = F),
    actionButton("btn_tab2_pca", "PCA", class = "btn btn-primary btn-block", disabled = F),
    actionButton("btn_tab2_umap", "UMAP", class = "btn btn-primary btn-block", disabled = F),
    actionButton("btn_tab2_clustering", "Clustering", class = "btn btn-primary btn-block", disabled = F),
    actionButton("btn_tab2_visualize", "Download Results", class = "btn btn-primary btn-block", disabled = F),
  ),
  div(
    id = "tab2_main_content",
    # Overview 内容
    div(
      id = "tab2_overview_panel",
      fluidRow(
        column(10,
               h4("Overview"),
               p("This section provides an overview of the multi-omics data analysis workflow."),
               img(src = "data/overview_img/tab2.overview.png", width = "800px", height = "800px"),
        ),
        column(2,
               img(src = "data/overview_img/logo.png", width = "200px", height = "100px")
        )
      )
    ),
  
    # Upload Data 内容，使用 fluidRow 实现左右布局
    div(
      id = "tab2_upload_panel",
      fluidRow(
        column(6,
               h4("Load example data"),
               actionButton("load_example", "Load Example and run", class = "btn btn-primary btn-hover"),
               
               h4("Upload your input data"),
               fileInput("user_data", "SingleCellExperiment Object (Accepted Format: .rds)"),
               
               h5("Analysis method: Seurat"),
               
               numericInput("min_feature", "Min nFeature:", value = 200),
               numericInput("min_cells", "Min nCells:", value = 3),
               
               textInput("Project_name", "Project Name:", value = "Project1"),
               actionButton("load_data", "Load Data and run", class = "btn btn-primary btn-hover"),
               actionButton("reset", "Reset", class = "btn btn-danger btn-hover")
        ),
        column(6,
               plotOutput("vln_plot_rna"),
               plotOutput("vln_plot_adt")
        )
      ),
      fluidRow(
        column(12,
               h4("RNA Data and ADT Data Information"),
               DTOutput("rna_data_table"),  # RNA Data 表格
               DTOutput("adt_data_table")   # ADT Data 表格
        )
      )
    ),
      
    # Normalize Data 内容
    div(
      id = "tab2_normalization_panel",
      h4("Data Preprocessing"),
      numericInput("num_var_feature", "Number of variable features:", value = 2000),
      selectInput("var_method", "Variable feature selection method:", choices = c("vst", "mean.var", "dispersion")),
      actionButton("run_variable_features", "Identify highly variable features", class = "btn btn-primary btn-hover"),
      br(),
      # output
      plotOutput("var_feature_plot", width=800),
    ),
      
    # PCA 内容
    div(
      id = "tab2_pca_panel",
      h4("Principal Component Analysis (PCA)"),
      actionButton("run_pca", "Run PCA", class = "btn-primary"),
      plotlyOutput("pca_plot", width=800),
      plotOutput("pca_viz", width=800),
    ),
  
    # UMAP 内容
    div(
      id = "tab2_umap_panel",
      h4("Uniform Manifold Approximation and Projection (UMAP)"),
      actionButton("run_umap", "Run UMAP", class = "btn-primary"),
      plotOutput("umap_plot_rna", width=800),
      plotOutput("umap_plot_adt", width=800),
      plotOutput("umap_plot_wnn", width=800)
    ),
    
    # Clustering 内容
    div(
      id = "tab2_clustering_panel",
      page_fillable(
        card(
          card_header("Clustering"),
          card_body(
            sliderInput("resolution_clustering", "Resolution:", value = 0.5, step = 0.1, min = 0.1, max = 3),
            actionButton("run_clustering", "Run Clustering", class = "btn-primary", style = "width: 20%"),
          ),
        ),
        card(
          card_header("RNA Clustering"),
          card_body(
            # 为RNA聚类图添加选择框
            selectInput("group_by_rna", "Group by for RNA:", choices = NULL),
            plotOutput("clustering_plot_rna", width=800),
          )
        ),
        card(
          card_header("ADT Clustering"),
          card_body(
            # 为ADT聚类图添加选择框
            selectInput("group_by_adt", "Group by for ADT:", choices = NULL),
            plotOutput("clustering_plot_adt", width=800),
          )
        ),
        card(
          card_header("WNN Clustering"),
          card_body(
            # 为WNN聚类图添加选择框
            selectInput("group_by_wnn", "Group by for WNN:", choices = NULL),
            plotOutput("clustering_plot_wnn", width=800),
          )
        )
      ),
      textOutput("clustering_result"),
    ),
    
    # Visualization 内容
    div(
      id = "tab2_visualize_panel",
      h4("Download Results"),
      downloadButton("download_rds", "Download", class = "btn btn-primary btn-hover"),
      # Metadata table display
      h4("Metadata Information"),
      tableOutput("metadata_table"),
      img(
        src = "data/overview_img/footer.png", 
        style = "position: fixed; bottom: 10px; right: 10px; width: 500px; height: 120px; z-index: 1000;"
      )
    )
  )
)