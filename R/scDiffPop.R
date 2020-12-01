

scDiffPop <- function(Sco, use.seurat.clusters = FALSE, find.markers = FALSE, find.pathways = FALSE) {
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

  results <- data.frame(group = rep(0,nrow(Tree)), enrichment=rep(0,nrow(Tree)), ES = rep(0, nrow(Tree)), p_adj = rep(0,nrow(Tree)), pvalue = rep(0,nrow(Tree)),
                        effect = rep(0,nrow(Tree)), lmpval = rep(0,nrow(Tree)), robust_stat = rep(0, nrow(Tree)),
                        robust_p = rep(0, nrow(Tree)), robust_fdr = rep(0, nrow(Tree)))
  GSEA_list <- list()
  responders.topgenes <- list()
  nonresponders.topgenes <- list()
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
    sco.sub <- Sco[,Sco$seurat_clusters %in% oldsubtree] %>%
      Seurat::FindVariableFeatures(selection.method = "vst", nfeatures = 2000)

    if(length(which(sco.sub@meta.data$binaryResponse == 0)) == 0
       || length(which(sco.sub@meta.data$binaryResponse == 1)) == 0) {
      gsea_result <- list()
      gsea_result$pvalue <- 1.0
      gsea_result$ES <- 0.0
      data[[i]]$GSEA <- gsea_result
      sco.sub <- Sco[,Sco$seurat_clusters %in% subtree]
      ncells <- nrow(sco.sub@meta.data)
      Counts[i,1] <- nrow(sco.sub@meta.data[sco.sub@meta.data$binaryResponse == 0,])
      Counts[i,2] <- nrow(sco.sub@meta.data[sco.sub@meta.data$binaryResponse == 1,])
      if(length(which(sco.sub@meta.data$binaryResponse == 0)) == 0) {
        results$enrichment[i] <- unique_phenotype[1]
        results$pvalue[i] <- 0
      }
      else{
        results$enrichment[i] <- unique_phenotype[2]
        results$pvalue[i] <- 0
      }
      next
    }

    ## Perform differential expression
    deg.curr = mydeg(sco.sub)
    #deg.curr = deg.curr[padj<0.1]

    ## Find markers in children
    seurat.clust.cell.1 = subtree
    Idents(sco.sub) = ifelse(sco.sub$seurat_clusters %in% seurat.clust.cell.1, 1 ,2) %>% as.factor
    markers.curr <- Seurat::FindMarkers(sco.sub, min.pct = 0.1, only.pos = TRUE, logfc.threshold = 0.25, ident.1 = 1)
    #markers.curr <- markers.curr[markers.curr$p_val < 0.05, ]
    markers.curr <- markers.curr[1:min(25, nrow(markers.curr)),]
    #markers.curr = FindAllMarkers(sco.sub, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25) ## this will give markers of both cell.1 (child1) and cell.2 (child2)
    #markers.top  = markers.curr %>% group_by(cluster) %>% top_n(n = 700, wt = avg_logFC)
    #term2gene  = markers.top %>%  as.data.table %>%
      #.[,.(term=cluster,gene=gene)] %>% as.data.frame
    print(nrow(markers.curr))
    term2gene <- matrix(0, nrow = nrow(markers.curr), ncol = 2) %>% as.data.frame
    term2gene[,1] <- 1; term2gene[,2] <- rownames(markers.curr)
    colnames(term2gene) <- c("term", "gene")

    ## GSEA

    if(nrow(markers.curr) == 0 || max(abs(deg.curr$log2FoldChange)) < 0.1) {
      gsea_result <- list()
      gsea_result$pvalue <- 1.0
      gsea_result$ES <- 0.0
      data[[i]]$GSEA <- gsea_result
      results$enrichment[i] <- unique_phenotype[1]
      results$pvalue[i] <- 1
    }
    else {
      #term2gene <- term2gene[term2gene$term == 1,]
      # GSEA
      geneList <- deg.curr$log2FoldChange
      names(geneList) <- deg.curr$gene
      geneList <- sort(geneList, decreasing = TRUE)
      print("MAX GENE: ")
      print(geneList[1])
      index <- which(rownames(Sco@assays$RNA@counts) == geneList[1])
      print(length(which(geneList > 0)))
      print("LIST SIZE:")
      print(length(geneList))
      print("MIN GENE:")
      print(geneList[length(geneList)])
      if(length(intersect(names(geneList), term2gene$gene)) == 0) {
        gsea_result <- list()
        gsea_result$pvalue <- 1.0
        gsea_result$ES <- 0.0
        data[[i]]$GSEA <- gsea_result
        results$enrichment[i] <- unique_phenotype[1]
        results$pvalue[i] <- 1
      }
      else {
        GSEA_list[[i]] <- clusterProfiler::GSEA(geneList, TERM2GENE = term2gene, pvalueCutoff = 1)

        gsea_result <- list()
        if(length(GSEA_list[[i]]@result$pvalue) == 0) {
          gsea_result$pvalue <- 1.0
          gsea_result$ES <- 0.0
          results$enrichment[i] <- unique_phenotype[1]
          results$pvalue[i] <- 1
        }
        else {
          gsea_result$ES <- GSEA_list[[i]]@result$enrichmentScore
          gsea_result$pvalue <- GSEA_list[[i]]@result$pvalue
          if(gsea_result$ES > 0) {
            results$enrichment[i] <- unique_phenotype[1]
            results$pvalue[i] <- gsea_result$pvalue
            results$ES[i] <- gsea_result$ES
          }
          else {
            results$enrichment[i] <- unique_phenotype[2]
            results$pvalue[i] <- gsea_result$pvalue
            results$ES[i] <- gsea_result$ES
          }
        }
        print(gsea_result)
        data[[i]]$GSEA <- gsea_result}
    }

    print(dim(sco.sub@assays$RNA@counts))

    genes.use <- which(rownames(sco.sub@assays$RNA@counts) %in% rownames(markers.curr))
    sco.sub@assays$RNA@counts <- sco.sub@assays$RNA@counts[genes.use,]
    sco.sub@assays$RNA@data <- sco.sub@assays$RNA@data[genes.use,]
    sco.sub <- FindVariableFeatures(sco.sub)

    print(dim(sco.sub@assays$RNA@counts))
    deg.curr <- mydeg(sco.sub)

    geneList <- deg.curr$log2FoldChange
    names(geneList) <- deg.curr$gene
    x <- markers.curr$avg_logFC
    names(x) <- rownames(markers.curr)
    y <- rep(0, length(x)); names(y) <- names(x)
    intsct1 <- which(names(geneList) %in% names(x))
    intsct2 <- which(names(x) %in% names(geneList))
    y[intsct2] <- geneList[intsct1]

    print(x)
    print(y)


    x[is.infinite(x)] <- 1000

    results$effect[i] <- 0
    results$lmpval[i] <- 1
    try({
    fit <- lm(y~x+0)
    print(summary(fit))
    fitsum <- summary(fit)
    results$effect[i] <- 25*fit$coefficients[1]
    print("EFFECT:")
    print(results$effect[i])
    results$lmpval[i] <- fitsum$coefficients[4]
    plot(x=x,y=y, xlab = "Marker l2FC", ylab=  "Phenotype l2FC")
    abline(lm(y~x+0), col = "red")})

    results$robust_stat[i] <- sum(x*y)
    null_dist <- permutation_test(sco.sub, 250, rownames(markers.curr), x)
    print(null_dist)
    pos_pval <- 1 - length(which(results$robust_stat[i] > null_dist))/250
    neg_pval <- 1 - length(which(results$robust_stat[i] < null_dist))/250
    results$robust_p[i] <- 2*min(pos_pval, neg_pval)

    hist(null_dist)
    abline(v = results$robust_stat[i], col = "red")

    sco.sub <- Sco[,Sco$seurat_clusters %in% subtree]
    sco.sub <- Seurat::FindVariableFeatures(sco.sub, selection.method = "vst", nfeatures = 2000)
    ncells <- nrow(sco.sub@meta.data)
    Counts[i,1] <- nrow(sco.sub@meta.data[sco.sub@meta.data$binaryResponse == 0,])
    Counts[i,2] <- nrow(sco.sub@meta.data[sco.sub@meta.data$binaryResponse == 1,])
    print(Counts[i,1])
    print(Counts[i,2])

    #Do differential expression
    try({
      deg.curr = mydeg(sco.sub)
      responders.topgenes[[i]] <- deg.curr[log2FoldChange > 0]$gene[1]
      nonresponders.topgenes[[i]] <- deg.curr[log2FoldChange < 0]$gene[1]})
  }

  results$p_adj <- p.adjust(results$pvalue, method = "fdr")
  results$robust_fdr <- p.adjust(results$robust_p, method = "fdr")

  p0 <- length(unique(Sco@meta.data$patient[Sco@meta.data$binaryResponse == 0]))
  p1 <- length(unique(Sco@meta.data$patient[Sco@meta.data$binaryResponse == 1]))


  xy <- layout_as_tree(G)
  V(G)$x <- xy[, 1]
  V(G)$y <- xy[, 2]

  pht1_size <- length(which(Sco@meta.data$binaryResponse == 0)); pht2_size <- length(which(Sco@meta.data$binaryResponse == 1))
  Counts <- rbind(c(pht1_size, pht2_size), Counts)

  l2total_counts <- log2(nrow(Sco@meta.data))
  V(G)$radius <- as.vector(log2(Counts[,1] + Counts[,2])/l2total_counts)
  print(V(G)$radius)
  print(Counts)

  Counts_c <- Counts
  Counts[,1] <- (Counts[,1]/p0) / (Counts[,1]/p0 + Counts[,2]/p1)
  Counts[,2] <- (Counts[,2]/p1) / (Counts_c[,1]/p0 + Counts[,2]/p1)

  V(G)$pht1 <- Counts[,1]
  V(G)$pht2 <- Counts[,2]

  print(V(G)$pht1); print(V(G)$pht2)

  #Save name for later
  name_clean <- V(G)$name
  V(G)$transparency <- rep(0.4, length(V(G)$name))

  for(i in 2:length(V(G)$name)) {
    s <- results$p_adj[i-1]
    if(s < 0.001) {
      V(G)$name[i] <- paste(V(G)$name[i], "***", sep="")
      V(G)$transparency[i] <- 1
    }
    else if(s < 0.01) {
      V(G)$name[i] <- paste(V(G)$name[i], "**", sep="")
      V(G)$transparency[i] <- 0.8
    }
    else if(s < 0.05) {
      V(G)$name[i] <- paste(V(G)$name[i], "*", sep="")
      V(G)$transparency[i] <- 0.6
    }
  }

  print(V(G)$transparency)

  #scaling constant
  c <- 0.2
  ncirc <- round(max(V(G)$radius)/min(V(G)$radius))
  V(G)$radius <- V(G)$radius * c

  print(igraph::as_data_frame(G, "vertices"))
  piechart_data <- igraph::as_data_frame(G, "vertices")
  visualizations <- list()

  p <- ggraph(G, "manual", x=  V(G)$x, y=V(G)$y)+ geom_edge_link()
  p <- p + geom_node_circle(aes(x0=x,y0=y,r=radius), colour = NA, show.legend = FALSE, data = piechart_data, fill="white")
  p <- p + geom_scatterpie(
    aes(x=x, y=y, r=radius, alpha = forcats::fct_inorder(name)),
    data = piechart_data,
    cols = c("pht1", "pht2"),
    colour = NA,
    legend_name = "Phenotype",
  ) + scale_alpha_manual(values = piechart_data$transparency, name = NULL, labels = NULL)
  p <- p + scale_fill_manual(values = c("turquoise", "hotpink1"), labels = c(unique_phenotype[1], unique_phenotype[2]))
  p <- p + geom_node_label(aes(label = name, angle = 90), repel = FALSE, nudge_y = 0.25, col = "midnightblue")
  p <- p + theme_graph()

  #out$counts <- Counts
  visualizations$pies <- p

  #mycols <- colorRampPalette(colors = c("blue", "grey", "red"))(100)
  #intensity <- round(50*results$effect + 50)
  #print(intensity)
  #intensity[intensity > 100] <- 100; intensity[intensity < 1] <- 1
  #piechart_data$intensity <- c(50, intensity)
  piechart_data$intensity <- c(1, abs(results$effect))
  piechart_data$intensity[piechart_data$intensity > 1] <- 1
  effect_col <- ifelse(results$effect < 0, "blue", "red")
  piechart_data$effect_col <- c("white", effect_col)

  V(G)$name <- name_clean
  p <- ggraph(G, "manual", x=  V(G)$x, y=V(G)$y) + geom_edge_link()
  p <- p + geom_node_circle(aes(x0=x,y0=y,r=radius), colour = "black", show.legend = FALSE, data = piechart_data, fill="white")
  p <- p + geom_node_circle(aes(x0=x,y0=y,r=radius, fill = forcats::fct_inorder(name), alpha = forcats::fct_inorder(name)), colour = NA, show.legend = FALSE, data = piechart_data) + scale_fill_manual(values = piechart_data$effect_col)
  p <- p + scale_alpha_manual(values = piechart_data$intensity)
  p <- p + geom_node_label(aes(label = name, angle = 90), repel = FALSE, nudge_y = 0.25, col = "midnightblue")
  p <- p + theme_graph()

  visualizations$effect <- p

  piechart_data$RS <- c(0, sign(results$robust_stat))
  piechart_data$RS <- ifelse(piechart_data$RS > 0, "turquoise", "hotpink1")
  piechart_data$RS[1] <- "white"
  piechart_data$RS_intensity <- rep(0.3, nrow(piechart_data))
  for(i in 2:length(name_clean)) {
    s <- results$robust_p[i-1]
    if(s < 0.01) {
      piechart_data$RS_intensity[i] <- 1
    }
    else if(s < 0.5) {
      piechart_data$RS_intensity[i] <- 0.75
    }
    else if(s < 0.1) {
      piechart_data$RS_intensity[i] <- 0.5
    }
  }

  V(G)$name <- name_clean
  p <- ggraph(G, "manual", x=  V(G)$x, y=V(G)$y) + geom_edge_link()
  p <- p + geom_node_circle(aes(x0=x,y0=y,r=radius), colour = "black", show.legend = FALSE, data = piechart_data, fill="white")
  p <- p + geom_node_circle(aes(x0=x,y0=y,r=radius, fill = forcats::fct_inorder(name), alpha = forcats::fct_inorder(name)), colour = NA, show.legend = FALSE, data = piechart_data) + scale_fill_manual(values = piechart_data$RS)
  p <- p + scale_alpha_manual(values = piechart_data$RS_intensity)
  p <- p + geom_node_label(aes(label = name, angle = 90), repel = FALSE, nudge_y = 0.25, col = "midnightblue")
  p <- p + theme_graph()

  visualizations$robust_stat <- p

  xy <- layout_as_tree(G)
  V(G)$x <- xy[, 1]
  V(G)$y <- xy[, 2]

  V(G)$resp <- c(" ", unlist(responders.topgenes))
  V(G)$nonresp <- c(" ", unlist(nonresponders.topgenes))

  p <- ggraph(G, "manual", x=  V(G)$x, y=V(G)$y) + geom_edge_link()
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

  results$lmpval <- p.adjust(results$lmpval, method = "fdr")
  results$group <- name_clean[-1]
  results <- results[order(results$p_adj),]
  out$results <- results

  return(out)
}


Mode <- function(x) {
  xu <- unique(x)
  return(xu[which.max(tabulate(match(x, xu)))])
}

SplitGroup <- function(Sco_sub, ixs) {
  if(length(ixs) == 2) {
    return(c(1,2))
  }

  Sco_sub <- Seurat::NormalizeData(Sco_sub)

  #Find variable features
  Sco_sub <- Seurat::FindVariableFeatures(Sco_sub, verbose = FALSE)

  #Scale data
  Sco_sub <- Seurat::ScaleData(Sco_sub)

  #Run PCA on the variable features. Get 50 dimensional embeddings
  Sco_sub <- Seurat::RunPCA(Sco_sub, verbose = FALSE)
  embeddings <- Seurat::Embeddings(object = Sco_sub, reduction = 'pca')[,1:50]

  pseudobulk <- matrix(0, nrow = 0, ncol = 50)
  for(i in ixs-1) {
    rxs <- which(Sco_sub@meta.data$seurat_clusters == i)
    pseudobulk <- rbind(pseudobulk, colMeans(embeddings[rxs,]))
  }
  print(pseudobulk)

  #Cluster via k means
  km <- kmeans(pseudobulk, centers = 2, iter.max = 10)$cluster

  return(km)
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
  print(dim(data.use))
  group.info[cells.1, "group"] <- "Group1"
  group.info[cells.2, "group"] <- "Group2"
  group.info[, "group"] <- factor(x = group.info[, "group"])
  group.info$wellKey <- rownames(x = group.info)
  dds1 <- DESeq2::DESeqDataSetFromMatrix(
    countData = data.use + 1, # Add 1 to the count
    colData = group.info,
    design = ~ group
  )
  dds1 <- DESeq2::DESeq(object = dds1, fitType = "local")
  res <- DESeq2::lfcShrink(dds = dds1, coef = 2)
  #res <- DESeq2::results(dds1)
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
  deseq.out = DESeq2DETest(data.use=exp.curr[,c(responders,nonresponders)], cells.1=responders, cells.2=nonresponders)
  deseq.dt = deseq.out %>%
    as.data.frame %>%
    mutate(gene=rownames(.)) %>%
    as.data.table() %>%
    .[order(pvalue)] %>%
    .[,padj:=p.adjust(pvalue, method="fdr")]
  deseq.dt
}


permutation_test <- function(Sco, iterations, markers, markers_avgFC) {
  null_dist <- rep(0, iterations)
  x <- markers_avgFC

  patients <- unique(Sco@meta.data$patient)
  patients_response <- sapply(patients, function(x) {
    row <- which(Sco@meta.data$patient == x)[1]
    return(Sco@meta.data$binaryResponse[row])
  })

  print(markers)

  print(patients_response)

  for(i in 1:iterations) {
    print("Rep: ")
    print(i)
    Sco.curr <- Sco
    patients.perm <- sample(1:length(patients), size = length(patients), replace = FALSE)
    for(j in 1:length(patients)) {
      Sco.curr@meta.data[Sco.curr@meta.data$patient == patients[j],]$binaryResponse <- patients_response[patients.perm[j] ]
    }

    print(unique(Sco.curr@meta.data$patient))
    print(sapply(unique(Sco.curr@meta.data$patient), function(x) {
      row <- which(Sco.curr@meta.data$patient == x)[1]
      return(Sco.curr@meta.data$binaryResponse[row])
    }))

    null_dist[i] <- 0
    try({
    deg.curr <- mydeg(Sco.curr)

    geneList <- deg.curr$log2FoldChange
    names(geneList) <- deg.curr$gene
    y <- rep(0, length(markers)); names(y) <- markers
    intsct1 <- which(names(geneList) %in% markers)
    intsct2 <- which(markers %in% names(geneList))
    y[intsct2] <- geneList[intsct1]

    x[is.infinite(x)] <- 1000
    null_dist[i] <- sum(x*y)
    print("VAL:")
    print(null_dist[i])
  })
  }

  null_dist[is.na(null_dist)] <- 0
  return(null_dist)
}


