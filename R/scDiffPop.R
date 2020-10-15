

scDiffPop <- function(Sco, use.seurat.clusters = TRUE) {
  if(!("cellType" %in% colnames(Sco@meta.data)) && !(use.seurat.clusters)) {
    stop("Seurat object meta data must have column 'cellType.")
  }
  if(!("response" %in% colnames(Sco@meta.data))) {
    stop("Seurat object meta data must have column 'response'.")
  }
  if(length(unique(Sco@meta.data$response)) != 2) {
    stop("Response column must have exactly two unique elements.")
  }

  if(use.seurat.clusters) {
    cell_types <- as.integer(Sco$seurat_clusters)
    Sco@meta.data$cellTypes <- sapply(cell_types, function(x) {paste("s", x, sep ="")})
  }

  cell_types <- unique(Sco@meta.data$cellType)
  Sco@meta.data$seurat_clusters <- sapply(Sco@meta.data$cellType, function(x) {
    return(which(cell_types == x) - 1)
  })

  group <- rep(1, length(cell_types))

  TreeMat <- matrix(1, nrow = length(cell_types), ncol = 1)
  counter <- 1

  while(length(unique(group)) != length(cell_types)) {
    #Get group with most clusters
    ixs <- which(group == Mode(group))
    print(group)
    cat("Current group: ", cell_types[ixs], "\n")
    split <- SplitGroup(Sco[,Sco$seurat_clusters %in% (ixs - 1)], ixs)

    if(length(unique(split)) > 1) {
      counter <- counter + 1
      subgroup <- group[ixs]
      subgroup[split == 2] <- counter
      counter <- counter + 1
      subgroup[split == 1] <- counter
      group[ixs] <- subgroup

      TreeMat <- cbind(TreeMat, group)
    }
    else {
      for(j in ixs) {
        counter <- counter + 1
        group[j] <- counter

        TreeMat <- cbind(TreeMat, group)
      }
    }
  }

  Tree <- matrix(nrow = 0, ncol = 2)

  counter <- 1

  for(i in 2:ncol(TreeMat)) {
    newvals <- unique(TreeMat[TreeMat[,i] > counter,i])
    for(j in newvals) {
      counter <- counter + 1
      ixs <- which(TreeMat[,i] == counter)
      newvec <- c(counter, TreeMat[ixs[1], i-1])
      Tree <- rbind(Tree, newvec)
    }
  }

  #Change names of leaves
  cntr <- 1
  for(ct in TreeMat[,ncol(TreeMat)]) {
    ix <- which(Tree[,1] == ct)
    Tree[ix,1] <- cell_types[cntr]
    cntr <- cntr+1
  }

  # Also give the columns a name
  colnames(Tree) <- c("Child", "Parent")
  rownames(Tree) <- c(1:nrow(Tree))

  print(Tree)

  TreeIG <- cbind(Tree[,2], Tree[,1])
  TreeIG <- as.matrix(TreeIG)

  G <- igraph::graph_from_edgelist(TreeIG)

  binaryResponse <- Sco$response
  unique_phenotype <- unique(Sco$response)
  binaryResponse[Sco$response == unique_phenotype[1]] <- 0
  binaryResponse[Sco$response == unique_phenotype[2]] <- 1
  binaryResponse <- as.integer(binaryResponse)
  Sco@meta.data$binaryResponse <- binaryResponse

  Counts <- matrix(0, nrow = nrow(Tree), ncol = 2)

  data <- list()
  data$subtree <- list()
  data$tree <- Tree

  # sco.sub = GSE145281[, GSE145281$seurat_clusters %in% c(10,11)]
  responders.enrichment = nonresponders.enrichment = list()
  responders.topgenes = nonresponders.topgenes = list()
  GSEA_list <- list()
  for (i in 1:nrow(Tree)) {
    data[[i]] <- list()
    cat("ITERATION: ", i, "\n")
    subtree <- as.vector(DFS(Tree, i, cell_types))
    data[[i]]$subtree <- subtree
    subtree <- which(cell_types %in% subtree) - 1
    print(subtree)
    old_ix <- as.integer(Tree[i,2]) - 1
    oldsubtree  <- which(cell_types %in% as.vector(DFS(Tree,old_ix,cell_types))) - 1
    print(oldsubtree)
    sco.sub = Sco[,Sco$seurat_clusters %in% oldsubtree] %>%
      Seurat::FindVariableFeatures(selection.method = "vst", nfeatures = 2000)

    if(length(which(sco.sub@meta.data$binaryResponse == 0)) == 0
       || length(which(sco.sub@meta.data$binaryResponse == 1)) == 0) {
      gsea_result <- list()
      gsea_result$pvalue <- 1.0
      gsea_result$ES <- 0.0
      data[[i]]$GSEA <- gsea_result
      next
    }
    ## Perform differential expression
    deg.curr = mydeg(sco.sub)
    #deg.curr = deg.curr[padj<0.1]

    ## Find markers in children
    seurat.clust.cell.1 = subtree
    Idents(sco.sub) = ifelse(sco.sub$seurat_clusters %in% seurat.clust.cell.1, 1 ,2) %>% as.factor
    markers.curr <- Seurat::FindMarkers(sco.sub, only.pos = TRUE, ident.1 = 1)
    markers.curr <- markers.curr[markers.curr$p_val < 0.05, ]
    #markers.curr = FindAllMarkers(sco.sub, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25) ## this will give markers of both cell.1 (child1) and cell.2 (child2)
    #markers.top  = markers.curr %>% group_by(cluster) %>% top_n(n = 700, wt = avg_logFC)
    #term2gene  = markers.top %>%  as.data.table %>%
      #.[,.(term=cluster,gene=gene)] %>% as.data.frame
    term2gene <- matrix(0, nrow = nrow(markers.curr), ncol = 2) %>% as.data.frame
    term2gene[,1] <- 1; term2gene[,2] <- rownames(markers.curr)
    colnames(term2gene) <- c("term", "gene")

    if(nrow(markers.curr) == 0) {
      gsea_result <- list()
      gsea_result$pvalue <- 1.0
      gsea_result$ES <- 0.0
      data[[i]]$GSEA <- gsea_result
      next
    }

    print("THROUGH MARKERS.CURR")

    # Make contingency table
    a <- length(which(deg.curr[log2FoldChange > 0]$gene %in% term2gene[term2gene[,1] == 1,2]))
    b <- length(which(deg.curr[log2FoldChange < 0]$gene %in% term2gene[term2gene[,1] == 1,2]))

    c <- length(deg.curr[log2FoldChange > 0]$gene) - a
    d <- length(deg.curr[log2FoldChange < 0]$gene) - b

    C <- matrix(c(a,b,c,d), nrow = 2, ncol = 2, byrow = TRUE)
    responders.enrichment[[i]] <- fisher.test(C, alternative = "greater")$p.value
    nonresponders.enrichment[[i]] <- fisher.test(C, alternative = "less")$p.value

    #term2gene <- term2gene[term2gene$term == 1,]
    # GSEA
    geneList <- deg.curr$log2FoldChange
    names(geneList) <- deg.curr$gene
    geneList <- sort(geneList, decreasing = TRUE)
    if(length(intersect(names(geneList), term2gene$gene)) == 0) {
      gsea_result <- list()
      gsea_result$pvalue <- 1.0
      gsea_result$ES <- 0.0
      data[[i]]$GSEA <- gsea_result
    }
    else {
    GSEA_list[[i]] <- clusterProfiler::GSEA(geneList, TERM2GENE = term2gene, pvalueCutoff = 1)

    gsea_result <- list()
    if(length(GSEA_list[[i]]@result$pvalue) == 0) {
      gsea_result$pvalue <- 1.0
      gsea_result$ES <- 0.0
    }
    else {
      gsea_result$ES <- GSEA_list[[i]]@result$enrichmentScore
      gsea_result$pvalue <- GSEA_list[[i]]@result$pvalue
    }
    print(gsea_result)
    data[[i]]$GSEA <- gsea_result}

    sco.sub = Sco[,Sco$seurat_clusters %in% subtree] %>%
      FindVariableFeatures(selection.method = "vst", nfeatures = 2000)
    ncells <- nrow(sco.sub@meta.data)
    Counts[i,1] <- nrow(sco.sub@meta.data[sco.sub@meta.data$binaryResponse == 0,])
    Counts[i,2] <- nrow(sco.sub@meta.data[sco.sub@meta.data$binaryResponse == 1,])
    print(Counts[i,1])
    print(Counts[i,2])

    #Do differential expression
    print("TRYING DIFFERENTIAL EXPRESSION")
    try({
      deg.curr = mydeg(sco.sub)
      responders.topgenes[[i]] <- deg.curr[log2FoldChange > 0]$gene[1]
      nonresponders.topgenes[[i]] <- deg.curr[log2FoldChange < 0]$gene[1]})
  }

  p0 <- length(unique(Sco@meta.data$patient[Sco@meta.data$binaryResponse == 0]))
  p1 <- length(unique(Sco@meta.data$patient[Sco@meta.data$binaryResponse == 1]))


  xy <- layout_as_tree(G)
  V(G)$x <- xy[, 1]
  V(G)$y <- xy[, 2]

  pht1_size <- length(which(Sco@meta.data$binaryResponse == 0)); pht2_size <- length(which(Sco@meta.data$binaryResponse == 1))
  Counts <- rbind(c(pht1_size, pht2_size), Counts)

  Counts_c <- Counts
  Counts[,1] <- (Counts[,1]/p0) / (Counts[,1]/p0 + Counts[,2]/p1)
  Counts[,2] <- (Counts[,2]/p1) / (Counts_c[,1]/p0 + Counts[,2]/p1)

  V(G)$pht1 <- Counts[,1]
  V(G)$pht2 <- Counts[,2]

  ix <- which(nonresponders.enrichment < 0.05 | responders.enrichment < 0.05)

  #Save name for later
  name_clean <- V(G)$name

  for(i in 2:length(V(G)$name)) {
    s <- data[[i-1]]$GSEA$pvalue
    if(s < 0.001) {
      V(G)$name[i] <- paste(V(G)$name[i], "***", sep="")
    }
    else if(s < 0.01) {
      V(G)$name[i] <- paste(V(G)$name[i], "**", sep="")
    }
    else if(s < 0.05) {
      V(G)$name[i] <- paste(V(G)$name[i], "*", sep="")
    }
  }

  visualizations <- list()

  p <- ggraph(G, "manual", x=  V(G)$x, y=V(G)$y) + geom_edge_link() #+ geom_node_point(aes(size = 2, col = names(effect)))
  #p <- p + geom_node_text(aes(x = x*1.005, y=y*1.005, label = name, angle = 90),
  p <- p +   geom_scatterpie(
    cols = c("pht1", "pht2"),
    data = igraph::as_data_frame(G, "vertices"),
    colour = NA,
    pie_scale = 0.75,
    legend_name = "Phenotype"
  )
  p <- p + scale_fill_manual(values = c("red", "blue"),
                             labels = c(unique_phenotype[1], unique_phenotype[2]))
  p <- p + geom_node_label(aes(label = name, angle = 90), repel = FALSE, nudge_y = 0.25, col = "midnightblue")
  p <- p + theme_graph()

  #out$counts <- Counts
  visualizations$pies <- p


  xy <- layout_as_tree(G)
  V(G)$x <- xy[, 1]
  V(G)$y <- xy[, 2]

  V(G)$resp <- c(" ", unlist(responders.topgenes))
  V(G)$nonresp <- c(" ", unlist(nonresponders.topgenes))

  p <- ggraph(G, "manual", x=  V(G)$x, y=V(G)$y) + geom_edge_link() #+ geom_node_point(aes(size = 2, col = names(effect)))
  #p <- p + geom_node_text(aes(x = x*1.005, y=y*1.005, label = name, angle = 90),
  #p <- p+ geom_node_label(aes(label = V(G)$Blood), repel = TRUE, col = "red")
  p <- p + geom_node_text(aes(label = V(G)$resp), repel = FALSE, nudge_x = -0.3, nudge_y = 0, col = "red")
  p <- p + geom_node_text(aes(label = V(G)$nonresp), repel = FALSE, nudge_x = -0.3, nudge_y = 0.15, col = "blue")
  p <- p + geom_node_point(size = 1)
  p <- p + geom_node_label(aes(label = name, angle = 90), repel = FALSE, nudge_y = 0.25, col = "midnightblue")
  #p <- p + geom_node_label(aes(label = name, angle = 90), repel = FALSE, nudge_y = 0.25, col = "midnightblue")
  p <- p + theme_graph()

  visualizations$genes <- p

  #Now do GSEA graph
  V(G)$name <- name_clean
  V(G)$ES <- sapply(data, function(x) {round(x$GSEA$ES,2)})
  V(G)$ES <- c("", V(G)$ES)
  p <- ggraph(G, "manual", x=  V(G)$x, y=V(G)$y) + geom_edge_link() #+ geom_node_point(aes(size = 2, col = names(effect)))
  #p <- p + geom_node_text(aes(x = x*1.005, y=y*1.005, label = name, angle = 90),
  #p <- p+ geom_node_label(aes(label = V(G)$Blood), repel = TRUE, col = "red")
  p <- p + geom_node_text(aes(label = V(G)$ES), repel = FALSE, nudge_x = -0.1, nudge_y = 0, col = "red")
  p <- p + geom_node_point(size = 1)
  p <- p + geom_node_label(aes(label = name, angle = 90), repel = FALSE, nudge_y = 0.25, col = "midnightblue")
  #p <- p + geom_node_label(aes(label = name, angle = 90), repel = FALSE, nudge_y = 0.25, col = "midnightblue")
  p <- p + theme_graph()

  visualizations$GSEA <- p

  out <- list()
  out$visualizations <- visualizations
  #data <- list()
  #data$GSEA <- GSEA_list; data$topgenes <- NULL
  #data$subtree <- subtree
  #data$phenotype <- data.frame(label = c(0,1), phenotype = c(unique_phenotype[1], unique_phenotype[2]))
  out$data <- data

  return(out)
}


Mode <- function(x) {
  xu <- unique(x)
  return(xu[which.max(tabulate(match(x, xu)))])
}

SplitGroup <- function(Sco_sub, ixs) {
  #Find variable features
  Sco_sub <- Seurat::FindVariableFeatures(Sco_sub, verbose = FALSE)

  #Scale data
  Sco_sub <- Seurat::ScaleData(Sco_sub)

  #Run PCA on the variable features. Get 50 dimensional embeddings
  Sco_sub <- Seurat::RunPCA(Sco_sub, verbose = FALSE)
  embeddings <- Seurat::Embeddings(object = Sco_sub, reduction = 'pca')[,1:50]

  #Cluster via k means
  km <- kmeans(embeddings, centers = 2, iter.max = 10)$cluster

  #Now we partition the clusters into one of two groups
  out <- sapply(ixs-1, function(x) {
    Mode(km[which(Sco_sub$seurat_clusters == x)])
  })

  return(out)
}

DFS <- function(Tree, node, cell_types) {
  if(node == 0) {
    return(cell_types)
  }

  leaves <- c()

  row <- Tree[node,]

  newparent <- row[1]

  allchild <- which(Tree[,2] == row[1])

  for(child in allchild) {
    leaves <- c(leaves, DFS(Tree, child))
  }

  if(length(allchild) == 0) {
    return(row[1])
  }

  return(leaves)
}

DESeq2DETest <- function(
  data.use,
  cells.1,
  cells.2,
  verbose = TRUE,
  ...
) {
  # if (!PackageCheck('DESeq2', error = FALSE)) {
  #   stop("Please install DESeq2 - learn more at https://bioconductor.org/packages/release/bioc/html/DESeq2.html")
  # }
  group.info <- data.frame(row.names = c(cells.1, cells.2))
  print(group.info)
  group.info[cells.1, "group"] <- "Group1"
  group.info[cells.2, "group"] <- "Group2"
  group.info[, "group"] <- factor(x = group.info[, "group"])
  group.info$wellKey <- rownames(x = group.info)
  dds1 <- DESeq2::DESeqDataSetFromMatrix(
    countData = data.use,
    colData = group.info,
    design = ~ group
  )
  dds1 <- DESeq2::DESeq(object = dds1, fitType = "local")
  res <- DESeq2::lfcShrink(dds = dds1, coef = 2)
  # to.return <- data.frame(p_val = res$pvalue, row.names = rownames(res))
  return(res)
}

mydeg <- function(sco.curr) {
  exp.curr1 = sco.curr@assays$RNA@counts[sco.curr@assays$RNA@var.features,]
  print(dim(exp.curr1))
  meta.dt1 = sco.curr@meta.data %>%
    as.data.table() %>%
    .[,.(binaryResponse=binaryResponse, patient=patient)]

  meta.curr = list()
  exp.curr2 = list()
  for(patient in unique(meta.dt1$patient)){
    inx = which(meta.dt1$patient==patient)
    exp.curr2[[patient]] = rowSums(as.matrix(exp.curr1[,inx]),na.rm=T)
    meta.curr[[patient]] = meta.dt1[inx[1],]
  }
  meta.dt = do.call(rbind, meta.curr)
  print(dim(meta.dt))
  exp.curr = t(do.call(rbind, exp.curr2))
  responders = meta.dt[binaryResponse==1]$patient
  nonresponders = meta.dt[binaryResponse==0]$patient
  print(responders)
  print(nonresponders)
  print(dim(exp.curr))
  print(dim(exp.curr[,c(responders,nonresponders)]))
  deseq.out = DESeq2DETest(data.use=exp.curr[,c(responders,nonresponders)], cells.1=responders, cells.2=nonresponders)
  deseq.dt = deseq.out %>%
    as.data.frame %>%
    mutate(gene=rownames(.)) %>%
    as.data.table() %>%
    .[order(pvalue)] %>%
    .[,padj:=p.adjust(pvalue, method="fdr")]
  deseq.dt
}


