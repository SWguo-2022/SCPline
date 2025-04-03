# 加载 Excel 文件
load_excel_file <- function(file_path, file_description) {
  incProgress(0.1, detail = paste("Loading", file_description, "..."))
  tryCatch({
    read_excel(file_path)
  }, error = function(e) {
    stop(paste("Error loading", file_description, ":", e$message))
  })
}

# 检查通道名称是否匹配
validate_channels <- function(panel_data, samp) {
  panel_channels <- trimws(panel_data$fcs_colname)
  fs_colnames <- trimws(colnames(samp[[1]]))
  
  channels_in_data <- panel_channels %in% fs_colnames
  
  if (!all(channels_in_data)) {
    missing_channels <- panel_channels[!channels_in_data]
    stop(paste("The following channels from the panel are not found in the FCS data:", paste(missing_channels, collapse = ", ")))
  }
}

# 主函数：加载用户数据
cyto_load_user_data <- function(mass_meta_data_path, mass_data_path, mass_fcs_zip_path, filetype) {
  # 加载元数据和面板数据
  meta_data <- load_excel_file(mass_meta_data_path, "Sample Metadata")
  panel_data <- load_excel_file(mass_data_path, "Panel Metadata")
  print("check file exists:")
  print(file.exists(mass_fcs_zip_path))
  temp_dir <- tempdir()
  print("list files of temp_dir:")
  list.files(temp_dir)
  # unlink(temp_dir, recursive = TRUE)  # 清理目录中的所有文件和子目录
  # dir.create(temp_dir)  # 重新创建临时目录
  zip::unzip(zipfile = mass_fcs_zip_path, exdir = temp_dir, overwrite = TRUE,junkpaths = F)
  
  fcs_files <- list.files(temp_dir, pattern = ".fcs$", full.names = TRUE, recursive=TRUE)
  print("fcs_files:")
  samp <- read.flowSet(files = fcs_files)
  
  # 验证通道名称匹配
  # validate_channels(panel_data, samp)
  
  incProgress(0.1, detail = "Preparing data...")
  
  # 预处理数据为 SCE 对象
  sce <- prepData(
    samp,
    panel = panel_data,
    md = meta_data,
    transform = TRUE,
    features = panel_data$fcs_colname,
    md_cols = list(
      file = "file_name", 
      id = "sample_id",
      group = "condition"
    )
  )
  
  return(list(sce = sce, meta_data = meta_data, panel_data = panel_data))
}

# 主函数：加载示例数据
LoadExampleMassData <- function(dataset = "./data/1.cytof_example_data", filetype = "fcs") {
  incProgress(0.1, detail = "Loading metadata information...")
  meta_data <- load_excel_file(file.path(dataset, "metadata", "sample_metadata.xlsx"), "Sample Metadata")
  panel_data <- load_excel_file(file.path(dataset, "metadata", "panel_metadata.xlsx"), "Panel Metadata")
  
  incProgress(0.1, detail = "Loading FCS data...")
  
  if (filetype == "fcs") {
    fcs_files <- list.files(file.path(dataset, "fcs_data"), pattern = "\\.fcs$", full.names = TRUE)
    
    if (length(fcs_files) == 0) {
      stop("No FCS files found in the specified directory.")
    }
    
    samp <- read.flowSet(files = fcs_files)
    
    # 验证通道名称匹配
    validate_channels(panel_data, samp)
    
    incProgress(0.1, detail = "Preparing data...")
    
    # 预处理数据为 SCE 对象
    sce <- prepData(
      samp,
      panel = panel_data,
      md = meta_data,
      transform = TRUE,
      features = panel_data$fcs_colname,
      md_cols = list(
        file = "file_name", 
        id = "sample_id",
        group = "condition"
      )
    )
    
    return(list(sce = sce, meta_data = meta_data, panel_data = panel_data))
  } else {
    stop("Unsupported filetype. Only 'fcs' is currently supported.")
  }
}

Auto_Gating <- function(sce, gating_strategy_file){
  # 检查 gating_strategy_file 是否存在
  if (!file.exists(gating_strategy_file)) {
    stop("The specified gating strategy file does not exist.")
  }
  
  # 加载门控策略
  gt <- gatingTemplate(gating_strategy_file)
  
  # 动态读取门控策略文件中的别名和父节点
  gating_strategy <- read.csv(gating_strategy_file, stringsAsFactors = FALSE)
  alias_names <- gating_strategy$alias
  parent_names <- gating_strategy$parent
  
  # 为每个细胞添加一个唯一的索引，用于追踪
  colData(sce)$cell_ind <- seq_len(ncol(sce))
  
  # 将 SCE 对象转换为 flowSet，按样本分割
  fs <- sce2fcs(sce, assay = "exprs", split_by = "sample_id", keep_cd = TRUE)
  
  # 构建 GatingSet 对象
  gs <- GatingSet(fs)
  
  # 应用门控策略
  gt_gating(gt, gs)
  
  # 提取门控频率，并存储为数据框
  gating_stats <- gs_pop_get_stats(gs, type = "percent")
  
  # 可选：打印门控统计信息
  print(as.data.frame(gating_stats))
  
  # 初始化一个列表来存储绘图对象
  plot_list <- list()
  
  # 获取 GatingSet 中所有群体的路径和别名映射
  pop_paths <- gs_get_pop_paths(gs, showHidden = TRUE)
  pop_aliases <- basename(pop_paths)
  alias_to_path <- setNames(pop_paths, pop_aliases)
  
  # 遍历门控别名，生成绘图并存储在列表中
  for (alias in alias_names)
  {
    # 检查别名是否存在于 GatingSet 中
    if (alias %in% pop_aliases) {
      # 获取对应的完整路径
      full_path <- alias_to_path[alias]
      
      # 使用 autoplot 生成门控图
      p <- autoplot(gs, full_path, bins = 1000, arrange = FALSE) +
        ggtitle(paste("Gating Plot for:", alias)) +
        axis_x_inverse_trans() + axis_y_inverse_trans()
      
      # 将绘图对象添加到列表中
      plot_list[[alias]] <- p
    } else {
      message("Alias '", alias, "' not found in GatingSet.")
    }
  }
  
  # 从 GatingSet 中提取最终门控后的数据
  # 动态获取最终门控别名
  final_alias <- alias_names[nrow(gating_strategy)]
  
  if (final_alias %in% pop_aliases) {
    final_path <- alias_to_path[final_alias]
  } else {
    stop("The final gating population '", final_alias, "' was not found in the GatingSet.")
  }
  
  # 获取通过最终门控的细胞数据
  fs_final <- gs_pop_get_data(gs, final_path)
  
  # 提取表达矩阵，并合并为单个数据框
  es_list <- lapply(fs_final, exprs)
  es <- do.call("rbind", es_list)
  
  # 使用 'cell_ind' 子集化 SCE 对象，保留通过门控的细胞
  sce_gated <- sce[, es[, "cell_ind"]]
  
  # 返回经过门控的 SCE 对象和绘图列表
  return(list(sce = sce_gated, plots = plot_list, stats = gating_stats))
}

tab1_server <- function(input, output, session) {
  tab1_current_page <- reactiveVal("tab1_overview_panel")
  observe({
    shinyjs::hide(selector = "#tab1_main_content > div")    
    shinyjs::show(tab1_current_page())
  })
  observeEvent(input$btn_cytof_overview, {
    tab1_current_page("tab1_overview_panel")
  })
  observeEvent(input$btn_cytof_upload, {
    tab1_current_page("tab1_upload_panel")
  })
  observeEvent(input$btn_cytof_visualization, {
    tab1_current_page("tab1_visualize_panel")
  })
  observeEvent(input$btn_cytof_normalization, {
    tab1_current_page("tab1_normalization_panel")
  })
  observeEvent(input$btn_cytof_auto_gating, {
    tab1_current_page("tab1_auto_gating_panel")
  })
  observeEvent(input$btn_cytof_analysis, {
    tab1_current_page("tab1_analysis_panel")
  })
  observeEvent(input$btn_cytof_download_results, {
    tab1_current_page("tab1_download_results_panel")
  })

  # 初始化流式数据对象
  meta_data_reactive <- reactiveVal(NULL)
  panel_data_reactive <- reactiveVal(NULL)
  sce_reactive <- reactiveVal(NULL)
  
  # 加载用户数据: 元数据、面板数据、FCS数据
  observeEvent(input$btn_load_user_data_meta, {
    # 验证文件输入
    req(input$mass_sample_data_path)
    req(input$mass_penal_data_path)
    req(input$mass_fcs_zip)
    
    # 获取上传文件路径
    sample_metadata_path <- input$mass_sample_data_path$datapath
    panel_metadata_path <- input$mass_penal_data_path$datapath
    fcs_zip_path <- input$mass_fcs_zip$datapath

    withProgress(message = "Loading user data...", value = 0, {
      incProgress(0.2, detail = "Processing files...")
      
      # 加载元数据和面板数据
      meta_data <- load_excel_file(sample_metadata_path, "Sample Metadata")
      panel_data <- load_excel_file(panel_metadata_path, "Panel Metadata")
      print("meta_data:")
      print(head(meta_data))
      print("panel_data:")
      print(head(panel_data))
      
      temp_dir <- tempdir()
      print("list files of temp_dir:")
      print(list.files(temp_dir))
      #
      
      # unlink(temp_dir, recursive = TRUE)  # 清理目录中的所有文件和子目录
      # dir.create(temp_dir)  # 重新创建临时目录
      print("list files of temp_dir:")
      print(list.files(temp_dir))
      zip::unzip(zipfile = fcs_zip_path, exdir = temp_dir, junkpaths = T)
      
      fcs_files <- list.files(temp_dir, pattern = ".fcs$", full.names = TRUE)
      print("fcs_files:")
      print(fcs_files)
      
      samp <- read.flowSet(files = fcs_files)
      
      # 验证通道名称匹配
      # validate_channels(panel_data, samp)
      
      incProgress(0.1, detail = "Preparing data...")
      
      # 预处理数据为 SCE 对象
      sce <- prepData(
        samp,
        panel = panel_data,
        md = meta_data,
        transform = TRUE,
        features = panel_data$fcs_colname,
        md_cols = list(
          file = "file_name", 
          id = "sample_id",
          group = "condition"
        )
      )

      # user_data <- cyto_load_user_data(sample_metadata_path, panel_metadata_path, fcs_zip_path, "fcs")
      
      # 如果加载失败，则停止处理
      # if (is.null(user_data)) {
      #   showNotification("Failed to load user data. Please check your files and try again.", type = "error")
      #   return()
      # }
      
      # 将加载的数据存储到反应式变量中
      sce_reactive(sce)
      meta_data_reactive(meta_data)
      panel_data_reactive(panel_data)
      
      # 显示成功通知
      incProgress(1, detail = "Done")
      showNotification("User data loaded successfully!", type = "message")
    })
    
    # 在分析步骤完成后调用
    showCompletionModal("Preprocess Data", "Data preprocessing completed!")
  })
  
  # 加载示例数据
  observeEvent(input$btn_load_example_data, {
    withProgress(message = "Loading example data...", detail = "This may take a few seconds.", value = 0, {
      example_data <- LoadExampleMassData()
      sce_reactive(example_data$sce)
      meta_data_reactive(example_data$meta_data)
      panel_data_reactive(example_data$panel_data)
    })
    showNotification("Example data loaded successfully!", type = "message")
    
    # 在分析步骤完成后调用
    showCompletionModal("Process Complete", "Example data loaded successfully!")
  })
  
  output$mass_data_path_text <- renderText({
    req(sce_reactive())
    sce_summary <- capture.output(print(sce_reactive()))
    HTML(paste("<pre>", paste(sce_summary, collapse = "\n"), "</pre>"))
  })
  
  output$meta_table <- renderTable({
    req(meta_data_reactive())
    meta_data_reactive()
  })
  
  output$panel_table <- renderTable({
    req(panel_data_reactive())
    panel_data_reactive()
  })
  
  # 1. Plot Cell Counts by Condition
  output$plot_counts <- renderPlot({
    req(sce_reactive())
    plotCounts(sce_reactive(), color_by = "condition", group_by = "sample_id", prop = FALSE)
  })
  
  # 2. Plot MDS
  output$plot_mds <- renderPlot({
    req(sce_reactive())
    pbMDS(sce_reactive(), color_by = "condition")
  })
  
  # 3. Expression Heatmap
  output$plot_heatmap <- renderPlot({
    req(sce_reactive())  # Ensure the SCE object is available
    plotExprHeatmap(sce_reactive(), bin_anno = TRUE, row_anno = TRUE)
  })
  
  # 4. Non-redundancy Score (NRS) Plot
  output$plot_nrs <- renderPlot({
    req(sce_reactive())
    plotNRS(sce_reactive(), features = type_markers(sce_reactive()), color_by = "condition", assay = "exprs")
  })
  
  # 5. Expression Boxplot
  output$plot_exprs <- renderPlot({
    req(sce_reactive())
    plotExprs(sce_reactive(), color_by = "condition")
  })

  # 6. Normalization Scatter Plot
  observeEvent(input$btn_normalization, {
    withProgress(message = "Performing normalization...", value = 0, {
      incProgress(0.5, detail = "normalization...")
      req(sce_reactive())
      # 6. Normalization Scatter Plot
      res <- normCytof(sce_reactive(), beads = "dvs", k = 50, assays = c("counts", "exprs"), overwrite = FALSE, plot = T)
      # 检查珠子和被移除事件的数量
      incProgress(0.1, detail = "Checking normalization results...")
      total_cells <- ncol(sce_reactive())
      num_beads <- ncol(res$beads)
      num_removed <- ncol(res$removed)
      incProgress(0.1, detail = "Building summary data frame...")
      # 构建摘要数据框
      summary_df <- data.frame(
        "Event Type" = c("Beads", "Removed"),
        "Count" = c(num_beads, num_removed),
        "Percentage (%)" = round(100 * c(num_beads / total_cells, num_removed / total_cells), 2)
      )
      incProgress(0.1, detail = "Building scatter plot...")
      output$summary_df <- renderTable(summary_df, rownames = FALSE)
      incProgress(0.1, detail = "Building scatter plot...")
      output$res_scatter <- renderPlot( res$scatter )
      incProgress(0.1, detail = "Building line plot...")
      # 7. Normalization Line Plot
      output$res_lines <- renderPlot(res$lines)
      
      sce <- res$data
      
      showNotification("Normalization completed!", type = "message")
    })
    
    # 在分析步骤完成后调用
    showCompletionModal("Process Complete", "Normalization completed successfully!")
  })
  
  # 8. auto gating
  # 定义函数：根据文件生成导航栏
  generate_nav_tabs <- function(file_path) {
    # 确保文件路径有效
    if (!file.exists(file_path)) {
      stop(paste("File does not exist:", file_path))
    }
    
    # 读取门控策略文件
    gating_strategy <- tryCatch(
      read.csv(file_path, stringsAsFactors = FALSE),
      error = function(e) {
        stop("Error reading gating strategy file:", e$message)
      }
    )
    
    # 检查 alias 列是否存在
    if (!"alias" %in% colnames(gating_strategy)) {
      stop("The gating strategy file does not contain an 'alias' column.")
    }
    
    # 获取所有别名，去重并删除无效值
    alias_names <- unique(na.omit(gating_strategy$alias))
    
    # 确保别名列表有效
    if (is.null(alias_names) || length(alias_names) == 0) {
      stop("Alias names are empty or invalid.")
    }
    
    # 动态生成 UI 导航栏
    output$dynamic_gating_tabs <- renderUI({
      panels <- lapply(alias_names, function(alias) {
        validate(
          need(alias != "", "Alias name cannot be empty")
        )
        nav_panel(
          title = paste("Auto Gating", alias),
          plotOutput(outputId = paste0("auto_gating_plot_", alias))
        )
      })
      
      # 检查生成的面板是否为空
      if (length(panels) == 0) {
        stop("No valid panels were generated for navigation.")
      }
      
      # 使用 do.call 解包 panels 列表
      do.call(navset_card_tab, panels)
    })
    
    # 返回别名列表，供后续绘图使用
    return(alias_names)
  }
  
  # 初始化导航栏（使用默认文件）
  observe({
    default_file <- "./data/1.cytof_example_data/metadata/gating_strategy.csv"
    if (file.exists(default_file)) {
      alias_names <- tryCatch(
        generate_nav_tabs(default_file),
        error = function(e) {
          showNotification(paste("Error initializing default navigation:", e$message), type = "error")
          return(NULL)
        }
      )
      
      if (!is.null(alias_names)) {
        lapply(alias_names, function(alias) {
          output[[paste0("auto_gating_plot_", alias)]] <- renderPlot({
            plot.new()
            text(0.5, 0.5, paste("No data available for", alias), cex = 1.2)
          })
        })
      }
    } else {
      showNotification("Default gating strategy file not found!", type = "error")
    }
  })
  
  # 监听文件上传事件
  observeEvent(input$gating_strategy_file, {
    req(input$gating_strategy_file)
    alias_names <- tryCatch(
      generate_nav_tabs(input$gating_strategy_file$datapath),
      error = function(e) {
        showNotification(paste("Error generating navigation for uploaded file:", e$message), type = "error")
        return(NULL)
      }
    )
    
    if (!is.null(alias_names)) {
      print(paste("Uploaded file aliases:", paste(alias_names, collapse = ", ")))
      lapply(alias_names, function(alias) {
        output[[paste0("auto_gating_plot_", alias)]] <- renderPlot({
          plot.new()
          text(0.5, 0.5, paste("Auto-gating in progress for", alias), cex = 1.2)
        })
      })
    }
  })
  
  # 监听按钮点击事件
  observeEvent(input$btn_auto_gating, {
    withProgress(message = "Performing auto-gating...", value = 0, {
      incProgress(0.1, detail = "Preparing data...")
      req(sce_reactive())  # 确保 SCE 对象存在
      
      # 判断文件来源：用户上传或默认文件
      gating_strategy_file <- if (is.null(input$gating_strategy_file)) {
        "./data/1.cytof_example_data/metadata/gating_strategy.csv"
      } else {
        input$gating_strategy_file$datapath
      }
      
      # 执行自动门控
      incProgress(0.5, detail = "Performing auto-gating...")
      gating_results <- Auto_Gating(sce_reactive(), gating_strategy_file)
      
      # 提取绘图列表
      plot_list <- gating_results$plots
      
      # 从门控策略文件动态获取别名
      gating_strategy <- read.csv(gating_strategy_file, stringsAsFactors = FALSE)
      alias_names <- gating_strategy$alias
      
      # 动态创建图表输出
      lapply(alias_names, function(alias) {
        output[[paste0("auto_gating_plot_", alias)]] <- renderPlot({
          req(plot_list[[alias]])
          plot_list[[alias]]
        })
      })
      
      # 显示门控统计信息
      stats_df <- as.data.frame(gating_results$stats)
      output$auto_gating_table <- renderTable(stats_df)
      
      # 写入门控结果到 CSV 文件
      output_file <- tempfile(fileext = ".csv")
      write.csv(stats_df, output_file, row.names = FALSE)
      
      # 提供下载功能
      output$download_gating_results <- downloadHandler(
        filename = function() {
          paste("gating_results_", Sys.Date(), ".csv", sep = "")
        },
        content = function(file) {
          file.copy(output_file, file)
        }
      )
      showNotification("Auto-gating completed and results available for download!", type = "message")
    })
    
    # 显示完成弹窗
    showCompletionModal("Process Complete", "Auto-gating completed successfully!")
  })
  # 9. Analysis
  observeEvent(input$btn_analysis_run, {
    withProgress(message = "Performing analysis...", value = 0, {
      incProgress(0.1, detail = "This step may take a few minutes...")
      req(sce_reactive())
      
      # 1. 聚类分析
      sce <- cluster(sce_reactive(), features = "type", xdim = 10, ydim = 10, maxK = 20, seed = 1234)
      message("Clustering done!")
      vi_cells <- 5000
      # 2. 降维
      message("Start running t-SNE and UMAP...")
      # sce <- runDR(sce, dr = "UMAP", cells = vi_cells, features = "type")
      sce <- runDR(sce, dr = "UMAP", features = "type")
      sce <- runDR(sce, dr = "TSNE", cells = vi_cells, features = "type")
      # 将处理后的 sce 存储到反应式值中
      sce_reactive(sce)
      
      # 看看UMAP降维数据存在
      print("看看UMAP降维数据是否存在:")
      print(head(reducedDim(sce, "UMAP"),40))
      print("Available reducedDims in sce_reactive():")
      print(names(reducedDims(sce_reactive())))
      
      # 更新分组选择的选项
      available_groups <- colnames(colData(sce)) # 获取所有列名, e.g. c("sample_id", "condition", "patient_id", "cluster_id")
      umap_groups = c("meta20", "meta10", available_groups)
      tsne_groups = c("meta20", "meta10", available_groups)
      updateSelectInput(session, "selected_umap_plot", choices = umap_groups, selected = umap_groups[1])
      updateSelectInput(session, "selected_tsne_plot", choices = tsne_groups, selected = tsne_groups[1])
      updateSelectInput(session, "selected_markers", choices = c("TSNE_markers", "UMAP_markers"))
      updateSelectInput(session, "selected_heatmap", choices = c("heatmap_meta10", "heatmap_meta20"))
      updateSelectInput(session, "selected_abundance_bar", choices = c("Abundance_barplot_meta10", "Abundance_barplot_meta20"))
      updateSelectInput(session, "selected_abundance_box", choices = c("Abundance_boxplot_meta10", "Abundance_boxplot_meta20"))
      updateSelectInput(session, "selected_expression", choices = c("Median_expr_meta10", "Median_expr_meta20"))
      
      # 在分析步骤完成后调用
      showCompletionModal("Process Complete", "Analysis completed successfully!")
    })
  })
  # 画图
  output$tsne_plot <- renderPlot({
    req(sce_reactive())
    req(input$selected_tsne_plot)
    plotDR(sce_reactive(), dr = "TSNE", color_by = input$selected_tsne_plot)
  })
  output$hhhumap_plot <- renderPlot({
    req(sce_reactive())
    req(input$selected_umap_plot)

    # Extract UMAP coordinates
    umap_coords <- reducedDim(sce_reactive(), "UMAP")
    
    # Print the dimensions and a preview of the UMAP coordinates
    print(dim(umap_coords))
    print(head(umap_coords))
    
    # Check for NA values
    if (any(is.na(umap_coords))) {
      showNotification("UMAP coordinates contain NA values", type = "error")
      return(NULL)  # Stop plotting if there are NA values
    }
    
    plotDR(sce_reactive(), dr = "UMAP", color_by = input$selected_umap_plot)
  })
  output$analysis_markers_plot <- renderPlot({
    req(input$selected_markers)
    req(sce_reactive())
    variable_markers <- rownames(sce_reactive())[order(matrixStats::rowVars(assay(sce_reactive(), "exprs")), decreasing = TRUE)][1:10]
    if (input$selected_markers == "TSNE_markers") {
      plotDR(sce_reactive(), dr = "TSNE", color_by = variable_markers)
    } else if (input$selected_markers == "UMAP_markers") {
      plotDR(sce_reactive(), dr = "UMAP", color_by = variable_markers)
    }
  })
  output$analysis_heatmap_plot <- renderPlot({
    req(input$selected_heatmap)
    req(sce_reactive())
    if (input$selected_heatmap == "heatmap_meta10") {
      k = "meta10"
    } else if (input$selected_heatmap == "heatmap_meta20") {
      k = "meta20"
    }
    variable_markers <- rownames(sce_reactive())[order(matrixStats::rowVars(assay(sce_reactive(), "exprs")), decreasing = TRUE)][1:10]
    plotExprHeatmap(sce_reactive(), k = k, features = variable_markers, by = "cluster_id", bars = TRUE, perc = TRUE)
  })
  output$analysis_abundance_barplot <- renderPlot({
    req(input$selected_abundance_bar)
    req(sce_reactive())
    if (input$selected_abundance_bar == "Abundance_barplot_meta10") {
      k = "meta10"
    } else if (input$selected_abundance_bar == "Abundance_barplot_meta20") {
      k = "meta20"
    }
    plotAbundances(sce_reactive(), k = k, by = "sample_id", group_by = "condition")
  })
  output$analysis_abundance_boxplot <- renderPlot({
    req(input$selected_abundance_box)
    req(sce_reactive())
    if (input$selected_abundance_box == "Abundance_boxplot_meta10") {
      k = "meta10"
    } else if (input$selected_abundance_box == "Abundance_boxplot_meta20") {
      k = "meta20"
    }
    plotAbundances(sce_reactive(), k = k, by = "cluster_id", group_by = "condition")
  })
  output$analysis_expression_plot <- renderPlot({
    req(input$selected_expression)
    req(sce_reactive())
    if (input$selected_expression == "Median_expr_meta10") {
      k = "meta10"
    } else if (input$selected_expression == "Median_expr_meta20") {
      k = "meta20"
    }
    variable_markers <- rownames(sce_reactive())[order(matrixStats::rowVars(assay(sce_reactive(), "exprs")), decreasing = TRUE)][1:10]
    plotPbExprs(sce_reactive(), k = k, features = variable_markers, assay = "exprs", group_by = "cluster_id")
  })
  
  # Download the processed data

  output$download_results <- downloadHandler(
    filename = function() {
      paste0("processed_data_", Sys.Date(), ".rds")
    },
    content = function(file) {
      req(sce_reactive())
      saveRDS(sce_reactive(), file = file)
    }
  )
  
}
