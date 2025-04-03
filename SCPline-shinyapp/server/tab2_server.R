tab2_server <- function(input, output, session) {
  # 默认显示 Overview 面板
  tab2_current_page <- reactiveVal("tab2_overview_panel")
  # 页面切换逻辑
  observe({
    # 隐藏所有页面
    shinyjs::hide(selector = "#tab2_main_content > div"    )    
    # 显示当前页面
    shinyjs::show(tab2_current_page())
  })
  observeEvent(input$btn_tab2_overview, {
    tab2_current_page("tab2_overview_panel")
  })
  observeEvent(input$btn_tab2_upload_data, {
    tab2_current_page("tab2_upload_panel")
  })
  observeEvent(input$btn_tab2_norm_feature, {
    tab2_current_page("tab2_normalization_panel")
  })
  observeEvent(input$btn_tab2_pca, {
    tab2_current_page("tab2_pca_panel")
  })
  observeEvent(input$btn_tab2_umap, {
    tab2_current_page("tab2_umap_panel")
  })
  observeEvent(input$btn_tab2_clustering, {
    tab2_current_page("tab2_clustering_panel")
  })
  observeEvent(input$btn_tab2_visualize, {
    tab2_current_page("tab2_visualize_panel")
  })
  
  # 处理上传的数据
  process_loaded_data <- function(srt_obj) {
    # 计算QC 指标
    srt_obj[['percent.mito']] <- PercentageFeatureSet(srt_obj, pattern = "^MT-")
    output$vln_plot_rna <- renderPlot({
      VlnPlot(srt_obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), pt.size = 0, ncol = 3)
    })
    output$vln_plot_adt <- renderPlot({
      VlnPlot(srt_obj, features = c("nFeature_ADT", "nCount_ADT"), pt.size = 0, ncol = 2)
    })
    
    # RNA data 信息表格
    output$rna_data_table <- renderDT({
      rna_data <- GetAssayData(srt_obj, assay = "RNA", layer = "counts")
      if (ncol(rna_data) > 10) {
        rna_data_subset <- as.data.frame(rna_data[1:10, 1:10])
      } else {
        rna_data_subset <- as.data.frame() # 空数据框
      }
      datatable(rna_data_subset, options = list(scrollX = TRUE))
    })
    
    # ADT data 信息表格
    output$adt_data_table <- renderDT({
      adt_data <- GetAssayData(srt_obj, assay = "ADT", layer = "counts")
      if (ncol(adt_data) > 10) {
        adt_data_subset <- as.data.frame(adt_data[1:10, 1:10])
      } else {
        adt_data_subset <- as.data.frame() # 空数据框
      }
      datatable(adt_data_subset, options = list(scrollX = TRUE))
    })
  }
  
  # 定义一个 reactiveVal 来存储 Seurat 对象
  seurat_obj <- reactiveVal()
  
  # 处理 load_example 按钮点击事件
  observeEvent(input$load_example, {
    withProgress(message = "Loading example data...", value = 0, {
      incProgress(0.5, detail = "Reading example data...")
      srt_obj <- readRDS("./data/sc_example_data/example_data.rds")
      if (class(srt_obj)[1] == "SingleCellExperiment") {
        srt_obj <- as.Seurat(srt_obj)
      } else if (class(srt_obj)[1] == "Seurat") {
        srt_obj
      } else {
        stop("Invalid data format. The RDS file should contain a SingleCellExperiment or Seurat object.")
      }
      
      # 删除包含 "res" 的列
      srt_obj@meta.data <- srt_obj@meta.data %>%
        select(-contains("res"))
      
      incProgress(0.5, detail = "Processing example data...")
      process_loaded_data(srt_obj)
      seurat_obj(srt_obj)
      
      # 在加载完成后调用
      showCompletionModal("Example Data Loaded", "Example data has been loaded successfully.")
    })
  })
  
  # 处理用户上传自己数据按钮点击事件
  observeEvent(input$user_data, {
    withProgress(message = "Loading user data...", value = 0, {
      incProgress(0.5, detail = "Reading user data...")
      srt_obj <- readRDS(input$user_data$datapath)
      if (class(srt_obj)[1] == "SingleCellExperiment") {
        srt_obj <- as.Seurat(srt_obj)
      } else if (class(srt_obj)[1] == "Seurat") {
        srt_obj
      } else {
        stop("Invalid data format. The RDS file should contain a SingleCellExperiment or Seurat object.")
      }
      incProgress(0.5, detail = "Processing user data...")
      process_loaded_data(srt_obj)
      seurat_obj(srt_obj)
      
      # 在加载完成后调用
      showCompletionModal("User Data Loaded", "User data has been loaded successfully.")
    })
  })
  
  # 处理run_variable_features 按钮点击事件
  observeEvent(input$run_variable_features, {
    withProgress(message = "Identifying variable features...", value = 0, {
      incProgress(0.5, detail = "Identifying variable features...")
      srt_obj <- req(seurat_obj())
      srt_obj <- NormalizeData(srt_obj)
      srt_obj <- FindVariableFeatures(srt_obj, selection.method = input$var_method, nfeatures = input$num_var_feature)
      top10_features <- head(VariableFeatures(srt_obj), 10)
      plot1 <- VariableFeaturePlot(srt_obj)
      plot2 <- LabelPoints(plot1, points = top10_features, repel = TRUE)
      output$var_feature_plot <- renderPlot(plot2)
      seurat_obj(srt_obj)
    })
    
    # 在分析步骤完成后调用
    showCompletionModal("Variable Features Complete", "Variable features have been identified successfully.")
  })
  
  # 处理run_pca 按钮点击事件
  observeEvent(input$run_pca, {
    withProgress(message = "Running PCA...", value = 0, {
      incProgress(0.5, detail = "Running PCA...")
      srt_obj <- req(seurat_obj())
      srt_obj <- ScaleData(srt_obj) %>% RunPCA(reduction.name = "pca_rna", npcs = 30)
      
      DefaultAssay(srt_obj) <- "ADT"
      adt_features <- rownames(srt_obj[['ADT']])
      srt_obj <- NormalizeData(srt_obj, normolization.method = 'CLR', margin = 2) %>% 
        ScaleData() %>%
        RunPCA(reduction.name = "pca_adt", features = adt_features, npcs = 30)
      
      output$pca_plot <- renderPlotly({
        # 确保 seurat_obj() 更新后再运行
        srt_obj <- req(seurat_obj())
        
        pca_data <- Embeddings(srt_obj, "pca_rna")
        pca_df <- data.frame(PC_1 = pca_data[, 1], PC_2 = pca_data[, 2])
        
        # 添加 Seurat 对象中的元数据，例如细胞类型
        meta_data <- srt_obj@meta.data
        pca_df$cluster <- meta_data$seurat_clusters  # 确保这个列名在你的meta数据中存在
        
        # 使用 Plotly 创建散点图
        plot_ly(pca_df, x = ~PC_1, y = ~PC_2, type = 'scatter', mode = 'markers',
                marker = list(color = ~cluster, colorscale = 'Viridis'),
                text = ~paste("Cluster:", cluster)) %>%  # 显示聚类信息
          layout(title = 'PCA Plot',
                 xaxis = list(title = 'PC 1'),
                 yaxis = list(title = 'PC 2'))
      })
      
      output$pca_viz <- renderPlot({
        # DimHeatmap(srt_obj, dims=1, cells=30, balanced=T, reduction="pca_rna")
        ElbowPlot(srt_obj, ndims = 20, reduction = "pca_rna")
      })
      
      seurat_obj(srt_obj)
    })
    # 在分析步骤完成后调用
    showCompletionModal("PCA Complete", "PCA analysis is finished and results are ready.")
  })
  
  # 处理 run_umap 按钮点击事件
  observeEvent(input$run_umap, {
    withProgress(message = "Running UMAP...", value = 0, {
      incProgress(0.2, detail = "prepare data for umap...")
      srt_obj <- req(seurat_obj())
      dims_rna <- length(srt_obj[["pca_rna"]]@stdev)  # 查看 pca_rna 中实际生成的 PCA 维度数
      dims_adt <- length(srt_obj[["pca_adt"]]@stdev)  # 查看 pca_adt 中实际生成的 PCA 维度数
      
      incProgress(0.1, detail = "Running UMAP for RNA data...")
      srt_obj <- RunUMAP(srt_obj, reduction = "pca_rna", dims = 1:dims_rna, reduction.name="umap_rna")
      incProgress(0.1, detail = "Running UMAP for adt data...")
      srt_obj <- RunUMAP(srt_obj, reduction = "pca_adt", dims = 1:dims_adt, reduction.name="umap_adt")
      pca_adt_dims <- dim(Embeddings(srt_obj, "pca_adt"))
      print(pca_adt_dims)
      print("finished umap_rna")
      
      incProgress(0.2, detail = "Running WNN UMAP...")
      srt_obj <- FindMultiModalNeighbors(srt_obj, reduction.list = list("pca_rna", "pca_adt"), 
                                         dims.list = list(1:dims_rna, 1:dims_adt),
                                         modality.weight.name = c("RNA.weight", "ADT.weight"))
      srt_obj <- RunUMAP(srt_obj, nn.name = "weighted.nn", reduction.name="wnnUMAP", reduction.key = "wnnUMAP_")
      print("finished wnnUMAP")
      
      seurat_obj(srt_obj)  # Update global reactive value
      
      # 仅展示 UMAP 图表，暂时不显示 groupby 选项
      incProgress(0.1, detail = "Rendering UMAP plots...")
      output$umap_plot_rna <- renderPlot({
        req(seurat_obj(), "umap_rna" %in% names(seurat_obj()@reductions))
        DimPlot(seurat_obj(), reduction = "umap_rna", label = TRUE) + ggtitle("RNA UMAP")
      })
      
      output$umap_plot_adt <- renderPlot({
        req(seurat_obj(), "umap_adt" %in% names(seurat_obj()@reductions))
        DimPlot(seurat_obj(), reduction = "umap_adt", label = TRUE) + ggtitle("ADT UMAP")
      })
      
      output$umap_plot_wnn <- renderPlot({
        req(seurat_obj(), "wnnUMAP" %in% names(seurat_obj()@reductions))
        DimPlot(seurat_obj(), reduction = "wnnUMAP", label = TRUE) + ggtitle("WNN UMAP")
      })
    })
      
    # 在分析步骤完成后调用
    showCompletionModal("Umap Complete", "UMAP analysis is finished and results are ready.")
  })
  
  # 处理 run_clustering 按钮点击事件
  observeEvent(input$run_clustering, {
    withProgress(message = "Running Clustering...", value = 0, {
      incProgress(0.2, detail = "Running Clustering...")

      srt_obj <- req(seurat_obj())
      dims_rna <- length(srt_obj[["pca_rna"]]@stdev)  # 查看 pca_rna 中实际生成的 PCA 维度数
      dims_adt <- length(srt_obj[["pca_adt"]]@stdev)  # 查看 pca_adt 中实际生成的 PCA 维度数
      
      srt_obj <- FindNeighbors(srt_obj, dims = 1:dims_rna, reduction = "pca_rna", graph.name = c("RNA_nn", "RNA_snn"))
      srt_obj <- FindNeighbors(srt_obj, dims = 1:dims_adt, reduction = "pca_adt", graph.name = c("ADT_nn", "ADT_snn"))
      
      srt_obj <- FindClusters(srt_obj, resolution = input$resolution_clustering, reduction = "umap_rna", graph.name = "RNA_snn")
      srt_obj <- FindClusters(srt_obj, resolution = input$resolution_clustering, reduction = "umap_adt", graph.name = "ADT_snn")
      print("finished clustering")
      seurat_obj(srt_obj)
      
      output$clustering_result <- renderText({
        paste("Clustering completed with", length(unique(Idents(srt_obj))), "clusters.")
      })
    })
    
    # 在分析步骤完成后调用
    showCompletionModal("Clustering Complete", "Clustering analysis is finished and results are ready.")
  })
  
  
  # 更新聚类图的选择框
  observe({
    req(seurat_obj())
    # 获取元数据中的所有列名
    meta_columns <- colnames(seurat_obj()@meta.data)
    # Find non-numeric columns
    non_numeric_columns <- meta_columns[!sapply(seurat_obj()@meta.data, is.numeric)]
    # 这里要为每个选择框更新选项,做一些筛选，比如去掉一些不需要的列，不该出现在ADT中的，不该出现在RNA中的
    group_by_rna_columns <- non_numeric_columns[!grepl("ADT|wnn", non_numeric_columns)]
    group_by_adt_columns <- non_numeric_columns[!grepl("RNA|wnn", non_numeric_columns)]
    group_by_wnn_columns <- non_numeric_columns[!grepl("RNA|ADT", non_numeric_columns)]
    
    # 更新每个选择框的选项
    updateSelectInput(session, "group_by_rna", choices = group_by_rna_columns)
    updateSelectInput(session, "group_by_adt", choices = group_by_adt_columns)
    updateSelectInput(session, "group_by_wnn", choices = group_by_wnn_columns)
  })
  
  # 为RNA聚类渲染图表
  output$clustering_plot_rna <- renderPlot({
    req(seurat_obj(), input$group_by_rna)
    DimPlot(seurat_obj(), reduction = "umap_rna", label = TRUE, group.by = input$group_by_rna, repel = TRUE) + ggtitle("RNA Clustering")
  })
  
  # 为ADT聚类渲染图表
  output$clustering_plot_adt <- renderPlot({
    req(seurat_obj(), input$group_by_adt)
    DimPlot(seurat_obj(), reduction = "umap_adt", label = TRUE, group.by = input$group_by_adt, repel = TRUE) + ggtitle("ADT Clustering")
  })
  
  # 为WNN聚类渲染图表
  output$clustering_plot_wnn <- renderPlot({
    # req(seurat_obj(), input$group_by_wnn)
    req(seurat_obj())
    print(head(seurat_obj()@meta.data))
    DimPlot(seurat_obj(), reduction = "wnnUMAP", label = TRUE, group.by = input$group_by_wnn, repel = TRUE) + ggtitle("WNN Clustering")
  })
  
  # 点击download_rds 时，下载 seurat 对象
  output$download_rds <- downloadHandler(
    filename = function() {
      "SingleCellExperiment.downloaded.rds"
    },
    content = function(file) {
      srt_obj <- req(seurat_obj())
      srt_obj <- as.SingleCellExperiment(srt_obj)
      saveRDS(srt_obj, file)
    }
  )
  
  output$metadata_table <- renderTable({
    req(seurat_obj())
    metadata <- seurat_obj()@meta.data
    if (ncol(metadata) > 10) {
      metadata_subset <- as.data.frame(metadata[1:10, 1:10])
    } else {
      metadata_subset <- as.data.frame() # 空数据框
    }
    metadata_subset
  })
}