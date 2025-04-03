tab1_ui <- page_sidebar(
  useShinyjs(),
  sidebar = sidebar(
    h4("Flow Cytometry"),
    actionButton("btn_cytof_overview", "Overview", class = "btn-primary"),
    actionButton("btn_cytof_upload", "Upload Data", class = "btn-primary"),
    actionButton("btn_cytof_normalization", "Normalization", class = "btn-primary"),
    actionButton("btn_cytof_auto_gating", "Auto Gating", class = "btn-primary"),
    actionButton("btn_cytof_visualization", "Statistic", class = "btn-primary"),
    actionButton("btn_cytof_analysis", "Main Analysis", class = "btn-primary"),
    actionButton("btn_cytof_download_results", "Download Results", class = "btn-primary"),
  ),
  div(
    id = "tab1_main_content",
    # tab1_overview_panel
    div(
      id = "tab1_overview_panel",
      fluidRow(
        column(10,
           id = "tab1_overview_panel",
           h4("Welcome to CyTOF workflow"),
           p("This is a shiny app for CyTOF data analysis."),
           p("Please click the buttons on the left to start."),
           img(src = "data/overview_img/tab1.overview.png", width = "800px", height = "800px"),
          ),
        column(2,
           img(src = "data/overview_img/logo.png", width = "200px", height = "100px")
        )
      )
    ),
    
    # tab1_upload_panel
    div(
      id = "tab1_upload_panel",
      fluidRow(
        column(12,
           card(
             card_header(h4("Load example data")),
             card_body(
               p("You can load example data to test the app."),
               actionButton("btn_load_example_data", "Run Example Data", class = "btn btn-primary btn-hover", style="width:20%")
             )
           )
        ),
        column(12,
           card(
             card_header(h4("Load User data")),
             card_body(
               fileInput("mass_sample_data_path", label = "Upload your sample metadata"),
               fileInput("mass_penal_data_path", label = "Upload your penal data"),
               fileInput("mass_fcs_zip", label = "Upload your FCS data ZIP file", accept = c(".zip")),
               actionButton("btn_load_user_data_meta", "Load and Run Data", class = "btn btn-primary btn-hover", style="width: 20%"),
             )
           )
        ),
        # 把load data的结果显示出来txt
        uiOutput("mass_data_path_text"),  # Change to uiOutput to accommodate the HTML format
        tableOutput("meta_table"),
        tableOutput("panel_table")
      )
    ),
    
    # tab1_normalization_panel
    div(
      id = "tab1_normalization_panel",
      fluidRow(
        column(12,
               card(
                 card_header(h4("Normalization")),
                 card_body(
                   p("This section provides normalization and feature selection."),
                   actionButton("btn_normalization", "Normalization", class = "btn btn-primary btn-hover", style="width:20%"),
                   tableOutput("summary_df")
                 )
               )
        ),
         navset_card_tab(
           nav_panel("Normalization Scatter Plot", plotOutput("res_scatter")),
           nav_panel("Normalization Line Plot", plotOutput("res_lines"))
         )
      )
    ),
    
    # tab1_visualize_panel
    div(
      id = "tab1_visualize_panel",
      # h4("Show some graphs"),
      navset_card_tab( 
        # Cell Counts Plot
        nav_panel("Cell Counts",  
                  div(style = "display: flex; justify-content: center;", 
                      plotOutput("plot_counts", height = "400px", width = "400px"))
        ),
        # MDS Plot
        nav_panel("MDS Plot", plotOutput("plot_mds")),
        
        # Expression Heatmap
        nav_panel("Expression Heatmap", plotOutput("plot_heatmap")),
        
        # Non-redundancy Score (NRS) Plot
        nav_panel("NRS Plot", plotOutput("plot_nrs")),
        
        # Exprs plot
        nav_panel("Exprs Plot", plotOutput("plot_exprs"))
      )
    ),
    
    # tab1_auto_gating_panel
    div(
      id = "tab1_auto_gating_panel",
      h4("Auto Gating"),
      p("Note: If you selected the example data in previous steps, you can proceed with auto-gating analysis by clicking the button below—no additional configuration needed.",
        "If you are using your own data, please upload your custom gating strategy file here to ensure the analysis aligns with your specific experimental setup."),
      
      fileInput("gating_strategy_file", label = "Upload your auto gating strategy file"),
      
      actionButton("btn_auto_gating", "Auto Gating", class = "btn btn-primary btn-hover"),
      
      div(style = "height: 150px; overflow-y: scroll;", 
          tableOutput("auto_gating_table")
      ),
      
      downloadButton("download_gating_results", "Download Gating Results"),
      
      # 动态创建导航栏和对应的绘图输出
      uiOutput("dynamic_gating_tabs")  # 将动态导航栏单独放置
    ),
    
    # tab1_analysis_panel
    div(
      id = "tab1_analysis_panel",
      h4("Analysis Results"),
      actionButton("btn_analysis_run", "Analysis Results", class = "btn btn-primary btn-hover"),
      br(),  # Adds a line break
      p(),
      navset_card_tab(
        nav_panel("TSNE Plot",
                  selectInput("selected_tsne_plot", "Select Plot:", choices = NULL),
                  plotOutput("tsne_plot")
        ),
        nav_panel("UMAP Plot", 
                  selectInput("selected_umap_plot", "Select Plot:", choices = NULL),
                  plotOutput("hhhumap_plot")
        ),
        nav_panel("Markers", 
                  selectInput("selected_markers", "Select Markers:", choices = c("TSNE_markers", "UMAP_markers")),
                  plotOutput("analysis_markers_plot")
        ),
        nav_panel("Heatmap", 
                  selectInput("selected_heatmap", "Select Heatmap:", choices = c("heatmap_meta20", "heatmap_meta10")),
                  plotOutput("analysis_heatmap_plot")
        ),
        nav_panel("Abundance Plot", 
                  selectInput("selected_abundance_bar", "Select Abundance:", choices = c("Abundance_barplot_meta20", "Abundance_barplot_meta10")),
                  plotOutput("analysis_abundance_barplot"),
                  selectInput("selected_abundance_box", "Select Abundance:", choices = c("Abundance_boxplot_meta20", "Abundance_boxplot_meta10")),
                  plotOutput("analysis_abundance_boxplot")
        ),
        nav_panel("Expression Plot", 
                  selectInput("selected_expression", "Select Expression:", choices = c("Median_expr_meta20", "Median_expr_meta10")),
                  plotOutput("analysis_expression_plot")
        )
      )
    ),
    
    # tab1_download_results_panel
    div(
      id = "tab1_download_results_panel",
      h4("Download Results"),
      downloadButton("download_results", "Download", class = "btn btn-primary btn-hover"),
      img(
        src = "data/overview_img/footer.png", 
        style = "position: fixed; bottom: 10px; right: 10px; width: 500px; height: 120px; z-index: 1000;"
      )
    )
  )
)
