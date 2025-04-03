library(QFeatures)
library(scp)
library(sva)

load_data <- function(peptide_file) {
  # 读取肽段数据
  peptide_data <- read.delim(peptide_file, sep = "\t", fileEncoding = "UTF-8")
  
  # 检查是否有 NA 列名
  if (anyNA(colnames(peptide_data))) {
    stop("Peptide data 文件中存在 NA 列名，请确保所有列名都是有效的。")
  }
  
  # 移除 peptide_data 中完全为空的列
  empty_cols <- sapply(peptide_data, function(x) all(is.na(x)))
  if (any(empty_cols)) {
    warning("Peptide data 中存在完全为空的列，这些列将被删除。")
    peptide_data <- peptide_data[, !empty_cols]
  }
  
  # 清理列名，确保其合法且唯一
  colnames(peptide_data) <- make.names(colnames(peptide_data), unique = TRUE)
  
  # 确定 ecol 的范围，基于列数
  ecol <- 3:ncol(peptide_data)
  
  # 使用 readQFeatures 读取数据
  cptac <- readQFeatures(peptide_data, quantCols = ecol, name = "peptides", fnames = "peptide")
  
  # 将 batch_label 和 celltype 数据合并或其他处理
  # 可根据需求处理 batch_label 和 celltype
  return(cptac)
}

load_example_data <- function() {
  withProgress(message = "Loading example data...", value = 0, {
    # 显示进度 25%
    incProgress(0.25, detail = "Reading peptide data...")
    
    peptide_file = "./data/3.mass_example_data/disrupted_peptide_level_data_pep_cell_0.5_6.txt"
    
    # 调用通用数据加载函数
    cptac <- load_data(peptide_file)
    
    # 显示进度 75%
    incProgress(0.75, detail = "Processing data...")
    
    # 显示完成进度 100%
    incProgress(1, detail = "Finalizing...")
    
    return(cptac)
  })
}

load_user_data <- function(input) {
  withProgress(message = "Loading user data...", value = 0, {
    # 显示进度 25%
    incProgress(0.25, detail = "Reading peptide data...")
    
    # 调用通用数据加载函数
    cptac <- load_data(peptide_file = input$mass_peptide_data_path$datapath)
    
    # 显示进度 75%
    incProgress(0.75, detail = "Processing data...")
    
    # 显示完成进度 100%
    incProgress(1, detail = "Finalizing...")
    return(cptac)
  })
}

process_data <- function(data, input) {
  selected_method <- input$mass_preprocess_method
  incProgress(1, detail = "Step 1: Replacing MaxQuant-encoded zeroes with NA values...")
  # Step 1: Replace MaxQuant-encoded zeroes with NA values
  cptac <- zeroIsNA(data, i = seq_along(data))
  print("in process_data function: ")

  # Step 2: Check and filter out data with more than 99% missing values
  incProgress(1, detail = "Step 2: Filtering out data with more than 99% missing values...")
  print("Step 2: Filtering out data with more than 99% missing values...")
  nNA(cptac, i = seq_along(cptac)) # Record results for later analysis if needed
  cptac <- filterNA(cptac, i = seq_along(cptac), pNA = 0.99)
  
  # Step 3: Count unique features at peptide and protein levels
  incProgress(1, detail = "Step 3: Counting unique features at peptide and protein levels...")
  print("Step 3: Counting unique features at peptide and protein levels...")
  if (!"peptide_counts" %in% colnames(colData(cptac))) {
    cptac <- countUniqueFeatures(cptac, i = "peptides", colDataName = "peptide_counts")
    cptac <- countUniqueFeatures(cptac, i = "peptides", groupBy = "protein", colDataName = "protein_counts")
  }
  
  # Step 4: Aggregate peptide data to protein level using multiple functions
  print("Step 4: Aggregating peptide data to protein level...")
  incProgress(1, detail = "Step 4: Aggregating peptide data to protein level...")
  # Feature aggregation with different functions
  if (selected_method == "sum") {
    cptac <- aggregateFeatures(cptac, i = "peptides", fcol = "protein", name = "proteins_colSums", fun = colSums)
    name = "proteins_colSums"
  } else if (selected_method == "mean") {
    cptac <- aggregateFeatures(cptac, i = "peptides", fcol = "protein", name = "proteins_colMeans", fun = colMeans)
    name = "proteins_colMeans"
  } else if (selected_method == "median") {
    cptac <- aggregateFeatures(cptac, i = "peptides", fcol = "protein", name = "proteins_colMedians", fun = colMedians)
    name = "proteins_colMedians"
  } else if (selected_method == "medianPolish") {
    cptac <- aggregateFeatures(cptac, i = "peptides", fcol = "protein", name = "proteins_medianPolish", fun = MsCoreUtils::medianPolish)
    name = "proteins_medianPolish"
  }
  # Step 5: Normalize columns and rows using median and mean centering
  incProgress(1, detail = "Step 5: Normalizing data by centering columns and rows...")
  print("Step 5: Normalizing data by centering columns and rows...")
  
  # ## Center columns with median ## Normalization
  # Step 6: Impute missing values using k-nearest neighbors (kNN)
  incProgress(1, detail = "Step 6: Imputing missing values using k-nearest neighbors (kNN)...")
  print("Step 6: Imputing missing values using k-nearest neighbors (kNN)...")
  
  # # 第[3]个 assay
  cptac <- sweep(cptac, i = name, MARGIN = 2, FUN = "-",
                 STATS = colMedians(assay(cptac[[name]]), na.rm = TRUE),
                 name = paste0(name, "_norm_col"))

  # Center rows with mean ## Normalization
  new_name = paste0(name, "_norm_col")
  new_name2 = paste0(name, "_norm")
  
  # 第[4]个 assay
  cptac <- sweep(cptac, i = new_name, MARGIN = 1, FUN = "-",
                 STATS = rowMeans(assay(cptac[[new_name]]), na.rm = TRUE),
                 name = paste0(name, "_norm"))

  assay(cptac[[new_name]]) <- assay(cptac[[new_name2]])

  # 把最后一个名字改成 proteins_colMedians_norm_imputed
  names(cptac)[4] <- paste0(names(cptac)[4], "_imputed")
  # names(cptac)[4] <- paste0(new_name2, "_imputed")
  # names(cptac)[5] <- "proteins_colMedians_norm"

  # Impute missing values using k-nearest neighbors (kNN)
  cptac <- impute(cptac, i = paste0(new_name2, "_imputed"),
                  method = "knn", k = 3, rowmax = 1, colmax = 1,
                  maxp = Inf, rng.seed = 1234)

  
  # No need to rename assays arbitrarily
  return(cptac)
}

run_pca <- function(cptac, selected_method) {
  name <- switch(selected_method,
                 "sum" = "proteins_colSums",
                 "mean" = "proteins_colMeans",
                 "median" = "proteins_colMedians",
                 "medianPolish" = "proteins_medianPolish")
  usename <- paste0(name, "_norm_imputed_batchC")
  if (!usename %in% names(cptac)) {
    usename <- paste0(name, "_norm_imputed")
  }
  
  # 提取数据矩阵并转换为稀疏矩阵
  library(Matrix)
  data_matrix <- assay(cptac[[usename]])
  data_matrix <- as(data_matrix, "dgCMatrix")
  
  if (anyNA(data_matrix) || any(is.infinite(data_matrix))) {
    stop("Data matrix contains NA or Inf values. Please clean your data.")
  }
  
  gene_variances <- apply(data_matrix, 1, var)
  if (all(gene_variances == 0)) {
    stop("All features have zero variance. PCA cannot be performed.")
  }
  
  gene_variances <- apply(data_matrix, 1, var)
  # print(summary(gene_variances))
  data_matrix <- data_matrix[gene_variances > 0, ]
  if (nrow(data_matrix) == 0) {
    stop("No features with non-zero variance found in the data matrix.")
  }
  
  # 创建 Seurat 对象
  seurat_obj <- CreateSeuratObject(counts = data_matrix, project = "Protein_Analysis", min.cells = 1, min.features = 10)
  seurat_obj[["protein"]] <- CreateAssayObject(data = data_matrix)
  DefaultAssay(seurat_obj) <- "protein"
  
  VariableFeatures(seurat_obj) <- rownames(data_matrix)
  seurat_obj <- ScaleData(seurat_obj, features = rownames(data_matrix), assay = "protein")
  seurat_obj <- RunPCA(seurat_obj, features = rownames(data_matrix), assay = "protein", verbose = FALSE)
  
  return(seurat_obj)
}

run_umap <- function(seurat_obj, cptac) {
  DefaultAssay(seurat_obj) <- "protein"
  seurat_obj <- FindNeighbors(seurat_obj, dims = 1:20)
  
  
  # **Specify the correct graph name in FindClusters**
  seurat_obj <- FindClusters(seurat_obj, resolution = 0.5, graph.name = "protein_snn")
  seurat_obj <- RunUMAP(seurat_obj, dims = 1:10)

  metadata_df <- as.data.frame(colData(cptac))

  if (nrow(metadata_df) != ncol(seurat_obj)) {
    stop("Error: The number of rows in metadata does not match the number of cells in Seurat object.")
  }
  
  seurat_obj <- AddMetaData(seurat_obj, metadata = metadata_df)
  return(seurat_obj)
}

tab3_server <- function(input, output, session) {
  current_page <- reactiveVal("tab3_overview3_panel")
  observe({
    shinyjs::hide(selector = "#tab3_main_content > div")
    shinyjs::show(current_page())
  })
  
  # 监听按钮点击事件，更新当前页面
  observeEvent(input$btn_mass_upload, {
    current_page("tab3_upload_panel")
  })
  observeEvent(input$btn_mass_preprocess, {
    current_page("tab3_preprocess_panel")
  })
  observeEvent(input$btn_mass_pca, {
    current_page("tab3_pca_panel")
  })
  observeEvent(input$btn_mass_umap, {
    current_page("tab3_umap_panel")
  })
  observeEvent(input$btn_mass_download, {
    current_page("tab3_download_panel")
  })
  observeEvent(input$btn_mass_overview3, {
    current_page("tab3_overview3_panel")
  })
  
  # 初始化
  cptac <- reactiveVal()
  metadata <- reactiveVal()
  seurat_obj <- reactiveVal()
  
  # 加载示例数据
  observeEvent(input$btn_mass_load_example_data, {
    example_data <- load_example_data()  # 加载示例数据
    if (inherits(example_data, "QFeatures")) {
      cptac(example_data)  # 确保 example_data 是 QFeatures 类型，然后保存到 cptac
      showNotification("Example data loaded successfully!", type = "message")
    } else {
      showNotification("Failed to load example data. Incorrect data type.", type = "error")
    }
    # 读取 metadata
    batch_label = read.csv("./data/3.mass_example_data/scope2_batchlabel.csv")
    celltype = read.csv("./data/3.mass_example_data/scope2_celltype.csv")
    metadata_data <- cbind(batch_label, celltype[, 2, drop = FALSE])
    
    # Update metadata reactiveVal
    metadata(metadata_data)
    # 更新 selectInput
    updateSelectInput(session, "mass_batch_effect_correction_field",
                      choices = c("batchlabel"),
                      selected = "batchlabel")
    # 在分析步骤完成后调用
    showCompletionModal("Load Example Data", "Example data loaded successfully!")
  })
  
  # 加载用户数据
  observeEvent(input$btn_mass_load_user_data, {
    user_data <- load_user_data(input)  # 加载用户数据
    if (inherits(user_data, "QFeatures")) {
      cptac(user_data)  
      showNotification("User data loaded successfully!", type = "message")
    } else {
      showNotification("Failed to load user data. Incorrect data type.", type = "error")
    }
    
    # 读取 metadata
    metadata_data <- read.csv(input$mass_metadata_path$datapath)
    
    # # Ensure sample names match
    # sample_names <- colnames(cptac())
    # if (!all(sample_names == metadata_data$SampleName)) {
    #   stop("Sample names in cptac and metadata_data do not match.")
    # }
    # rownames(metadata_data) <- metadata_data$SampleName
    # metadata_data$SampleName <- NULL
    # 
    # # Add metadata to cptac
    # colData(cptac()) <- cbind(colData(cptac()), metadata_data)
    
    metadata(metadata_data)
    # 更新 selectInput
    updateSelectInput(session, "mass_batch_effect_correction_field",
                      choices = colnames(metadata_data),
                      selected = colnames(metadata_data)[1])
    # 在分析步骤完成后调用
    showCompletionModal("Load User Data", "User data loaded successfully!")
  })
  
  # 预处理数据
  observeEvent(input$btn_mass_preprocess_data, {
    req(cptac())  # 验证 cptac 是否存在并有效
 
    withProgress(message = "Preprocessing data...", value = 0, {
      showNotification("Starting data preprocessing...", type = "message")
      
      # 处理数据，传递 reactiveVal 中的 cptac 对象
      processed_data <- process_data(cptac(), input)  # 将 cptac 中的数据传递给 process_data 处理
      
      metadata_data <- metadata()
      # Add metadata to 
      colData(processed_data) <- cbind(colData(processed_data), metadata_data)

      # 将预处理后的数据保存回 cptac
      cptac(processed_data)  
      
      incProgress(1, detail = "Data preprocessing completed")
      showNotification("Data preprocessing completed!", type = "message")
    })
    
    # 在分析步骤完成后调用
    showCompletionModal("Preprocess Data", "Data preprocessing completed!")
  })
  
  # 监控correct batch effect按钮
  observeEvent(input$btn_mass_batch_effect_correction, {
    req(cptac())  
    req(metadata())
    req(input$mass_batch_effect_correction_field)
    
    cptac <- cptac()
    selected_method <- input$mass_preprocess_method

    name <- switch(selected_method,
                   "sum" = "proteins_colSums",
                   "mean" = "proteins_colMeans",
                   "median" = "proteins_colMedians",
                   "medianPolish" = "proteins_medianPolish")

    withProgress(message = "Batch Effect Correction...", value = 0, {
      showNotification("Starting batch effect correction...", type = "message")
      
      # 提取用户选择的字段
      selected_field <- input$mass_batch_effect_correction_field
      
      # 从 metadata 中提取 batch_label 和 celltype 数据
      meta_data <- metadata()
      if (!selected_field %in% colnames(meta_data)) {
        showNotification("Selected field not found in metadata.", type = "error")
        stop("Invalid field selected for batch effect correction.")
      }
      
      # 确保元数据与 colData(cptac) 行数一致
      if (nrow(meta_data) != nrow(colData(cptac()))) {
        showNotification("Metadata row count does not match colData(cptac).", type = "error")
        stop("Metadata and colData(cptac) row mismatch.")
      }
      
      # print(colnames(cptac))
      # Step 7: Add batch and celltype data for batch effect correction
      incProgress(1, detail = "Step 7: Adding batch and cell type data for batch effect correction...")
      print("Step 7: Adding batch and cell type data for batch effect correction...")
      
      # 动态添加 batch_label 和 celltype 到 colData(cptac)
      colData(cptac) <- cbind(
        colData(cptac),
        batchlabel = meta_data[, selected_field, drop = T],
        celltype = meta_data[, "celltype", drop = T]
      )
      
      # print("8.after add batch and celltype:")
      # print(colData(cptac))
      # Step 8: Correct batch effects using the ComBat function
      incProgress(1, detail = "Step 8: Correcting batch effects using the ComBat method...")
      # print("Step 8: Correcting batch effects using the ComBat method...")
      cptac_sub <- getWithColData(cptac, paste0(name, "_norm_imputed"))
      
      batch <- colData(cptac_sub)$batchlabel
      model <- model.matrix(~ celltype, data = colData(cptac_sub))
      
      assay(cptac_sub) <- ComBat(dat = assay(cptac_sub),
                                 batch = batch,
                                 mod = model)
      # print("Combat...")
      # print(cptac)
      # Step 9: Link the corrected data back to the main cptac object
      incProgress(1, detail = "Step 9: Linking corrected data back to the main object...")
      # print("Step 9: Linking corrected data back to the main object...")
      cptac <- addAssay(cptac, y = cptac_sub, name = paste0(name, "_norm_imputed_batchC"))
      # Link the corrected data back to the main cptac object
      cptac <- addAssayLinkOneToOne(cptac, from = paste0(name,"_norm_imputed"),
                                    to = paste0(name,"_norm_imputed_batchC"))
      # print("finished...")
      # print(cptac)
      # Update the reactiveVal with the corrected data
      cptac(cptac)
      
      incProgress(1, detail = "Batch effect correction completed")
      showNotification("Batch effect correction completed!", type = "message")
    })
    
    # 在分析步骤完成后调用
    showCompletionModal("Batch Effect Correction", "Batch effect correction completed!")
  })
  
  # PCA
  observeEvent(input$btn_mass_pca_button, {
    req(cptac())
    req(input$mass_preprocess_method)
    
    withProgress(message = "PCA...", value = 0, {
      showNotification("Starting PCA analysis...", type = "message")
      
      # PCA 分析
      seurat <- run_pca(cptac(), input$mass_preprocess_method)
      
      # 获取 PCA 坐标并更新图表
      pca_embeddings <- Embeddings(seurat, "pca")
      mass_pca_data <- as.data.frame(pca_embeddings[, 1:2])
      
      colnames(mass_pca_data) <- c("PC_1", "PC_2")
      mass_pca_data$Sample <- rownames(mass_pca_data)
      
      output$mass_pca_plot <- renderPlotly({
        plot_ly(
          data = mass_pca_data,
          x = ~PC_1,
          y = ~PC_2,
          text = ~Sample,
          type = 'scatter',
          mode = 'markers'
        )
      })
      seurat_obj(seurat)  # 保存到 seurat_obj 以供 UMAP 使用
      
      incProgress(1, detail = "PCA completed")
      showNotification("PCA completed!", type = "message")
    })
    
    # 在分析步骤完成后调用
    showCompletionModal("PCA Analysis", "PCA analysis completed!")
  })
  
  # UMAP 
  observeEvent(input$btn_mass_umap_button, {
    req(seurat_obj())  # 确认 PCA 已运行并生成 seurat_obj
    req(cptac())  # 确认 cptac 已加载
    
    withProgress(message = "UMAP...", value = 0, {
      showNotification("Starting UMAP analysis...", type = "message")
      
      # 使用保存的 seurat_obj 进行 UMAP
      seurat <- run_umap(seurat_obj(), cptac())
      seurat_obj(seurat)  # 更新 seurat_obj
      
      incProgress(1, detail = "UMAP analysis completed")
      showNotification("UMAP analysis completed!", type = "message")
      # print("UMAP analysis completed!")
    })
    # Update the selectInput choices for 'mass_umap_group_by'
    observe({
      req(seurat_obj())
      metadata_columns <- colnames(seurat_obj()@meta.data)
      metadata_columns <- metadata_columns[!metadata_columns %in% c("nCount_RNA", "nFeature_RNA","protein_counts","peptide_counts","X")]
      updateSelectInput(session, "mass_umap_group_by",
                        choices = metadata_columns,
                        selected = metadata_columns[1])
    })
    # 在分析步骤完成后调用
    showCompletionModal("UMAP Analysis", "UMAP analysis completed!")
  })
  
  
  # 响应用户更改 group_by 选项并实时更新 UMAP 图
  observe({
    req(seurat_obj())
    req(input$mass_umap_group_by)
    
    group_by <- input$mass_umap_group_by
    
    output$umap_plot <- renderPlot({
      # 确保 UMAP 已经运行
      if ("umap" %in% names(seurat_obj()@reductions)) {
        DimPlot(seurat_obj(), reduction = "umap", group.by = group_by, label = TRUE, pt.size = 2)
      } else {
        text(0, 0, "UMAP not yet available. Please run the UMAP analysis.")
      }
    })
  })
  
  # 保存数据
  output$btn_mass_download_results <- downloadHandler(
    filename = function() {
      "MASS_SingleCellExperiment.Object.rds"
    },
    content = function(file) {
      srt_obj <- req(seurat_obj())
      srt_obj <- as.SingleCellExperiment(srt_obj)
      saveRDS(srt_obj, file)
    }
  )
}
