ancova.test <- function(df, factor, covar, trim = TRUE) {
  ##Prepare a blank data.frame for results
  result <-
    data.frame(
      "ID" = colnames(df),
      "p.value" = NA,
      "Adj.p.value" = NA,
      stringsAsFactors = FALSE
    )
  for (i in levels(factor)) {
    result[i] = NA
  }
  ##Prepare the formula used by ANCOVA
  fml_w_factor = as.formula(paste0("df[, i] ~ factor + ", paste(covar, collapse = " + ")))
  fml_wo_factor = as.formula(paste0("df[, i] ~ ", paste(covar, collapse = " + ")))
  ##Perform ANCOVA
  for (i in 1:ncol(df)) {
    if (sum(df[, i], na.rm = T) == 0) {
      next()
    }
    temp <-
      anova(lm(fml_w_factor), lm(fml_wo_factor))
    result[i, 2] <- temp$Pr[2]
    for (j in 4:ncol(result)) {
      result[i, j] = mean(subset(df[, i], subset = factor == levels(factor)[j - 3]))
    }
  }
  if (trim) {
    result <- result[complete.cases(result[, 2]), ]
  }
  result$Adj.p.value <- p.adjust(result$p.value, method = "fdr")
  
  return(result)
}

anova.test <- function(df, factor, trim = TRUE) {
  ##Prepare a blank data.frame for results
  result <-
    data.frame(
      "ID" = colnames(df),
      "p.value" = NA,
      "Adj.p.value" = NA,
      stringsAsFactors = FALSE
    )
  for (i in levels(factor)) {
    result[i] = NA
  }
  ##Perform ANOVA
  for (i in 1:ncol(df)) {
    if (sum(df[, i], na.rm = T) == 0) {
      next()
    }
    temp <- anova(lm(df[, i] ~ factor))
    result[i, 2] <- temp$Pr[1]
    for (j in 4:ncol(result)) {
      result[i, j] = mean(subset(df[, i], subset = factor == levels(factor)[j - 3]))
    }
  }
  if (trim) {
    result <- result[complete.cases(result[, 2]), ]
  }
  result$Adj.p.value <- p.adjust(result$p.value, method = "fdr")
  
  return(result)
}

b.diver <-
  function(df,
           factor,
           method = "bray",
           pal = "Dark2",
           lab.x = NULL,
           lab.legend = NULL)
  {
    library(vegan)
    library(ggplot2)
    library(ggpubr)
    library(foreach)
    ##Split the input data.frame by the factor
    dist_list <-
      foreach(i = levels(factor), .combine = "rbind") %do% {
        df_sub <- subset(df, subset = factor == i)
        dist_sub <-
          as.data.frame(as.numeric(vegdist(df_sub, method = method)))
        dist_sub["group"] = i
        dist_sub
      }
    colnames(dist_list) <- c("dist", "group")
    
    result <- list()
    
    if (length(levels(factor)) == 2) {
      ##While there are two groups in the factor
      result["test"] <-
        list(wilcox.test(dist ~ group, data = dist_list))
      test <-
        data.frame(
          group1 = levels(factor)[1],
          group2 = levels(factor)[2],
          p = result$test$p.value,
          y.position = max(dist_list$dist) * 1.02,
          p.signif = NA
        )
      if (test$p >= 0.05) {
        test$p.signif = "ns"
      } else if (test$p >= 0.01) {
        test$p.signif = "*"
      } else if (test$p >= 0.001) {
        test$p.signif = "**"
      } else {
        test$p.signif = "***"
      }
      
      result["plot"] <- list(
        ggplot(dist_list, aes(x = group, y = dist)) +
          scale_fill_brewer(palette = pal) +
          geom_boxplot(aes(fill = group)) +
          labs(x = lab.x, y = "Beta Diversity", fill = lab.legend) +
          stat_pvalue_manual(test, label = "p.signif", hide.ns = T)
      )
    } else {
      ##While there are three or more groups in the factor
      result["test"] <-
        list(kruskal.test(dist ~ group, data = dist_list))
      
      result["plot"] <- list(
        ggplot(dist_list, aes(x = group, y = dist)) +
          scale_fill_brewer(palette = pal) +
          geom_boxplot(aes(fill = group)) +
          labs(x = lab.x, y = "Beta Diversity", fill = lab.legend)
      )
    }
    
    return(result)
  }

build.network <- function(cor_matrix, p_matrix) {
  require("igraph")
  require("reshape2")
  
  temp <- cor_matrix
  temp[lower.tri(temp, diag = T)] = NA
  cor <- melt(temp)
  name <- rep(colnames(cor_matrix), ncol(cor_matrix))
  temp <- p_matrix
  temp[lower.tri(temp, diag = T)] = NA
  p <- melt(temp)
  name <- rep(colnames(p_matrix), ncol(p_matrix))
  cor_p_table <- cbind(name, cor, p)
  cor_p_table <- cor_p_table[, -2]
  colnames(cor_p_table) <- c("ID1", "cor", "ID2", "p")
  cor_p_table <- cor_p_table[, c("ID1", "ID2", "cor", "p")]
  cor_p_table$p <- p.adjust(cor_p_table$p, method = "BH")
  cor_p_table <- subset(cor_p_table, subset = cor_p_table$p < 0.05)
  cor_p_table <-
    subset(cor_p_table, subset = abs(cor_p_table$cor) < 0.9999)
  cor_p_table <-
    subset(cor_p_table, subset = abs(cor_p_table$cor) > 0.6)
  network_graph <-
    graph_from_edgelist(as.matrix(cor_p_table[, 1:2]), directed = F)
  E(network_graph)$weight = abs(cor_p_table$cor)
  network_graph <- simplify(network_graph)
  
  fine = 100
  palette = colorRampPalette(c("#56B4E9", "#E69F00"))
  
  eig <- eigen_centrality(network_graph)$vector
  eigColor <- palette(fine)[as.numeric(cut(eig, breaks = fine))]
  bet <- betweenness(network_graph)
  betColor <- palette(fine)[as.numeric(cut(bet, breaks = fine))]
  
  result <- list()
  result["networkGraph"] <- list(network_graph)
  result["eig"] <- list(eig)
  result["eigColor"] <- list(eigColor)
  result["bet"] <- list(bet)
  result["betColor"] <- list(betColor)
  
  return(result)
}

dataset.subset <- function (dataset, factor, subset) {
  library(foreach)
  
  result <- list()
  
  for (i in names(dataset)) {
    if (is.list(dataset[[i]])) {
      result[[i]] <- foreach (j = dataset[[i]]) %do% {
        subset(j, subset = factor %in% subset)
      }
      names(result[[i]]) = names(dataset[[i]])
    } else {
      result[i] <- list(subset(dataset[i], subset = factor %in% subset))
    }
  }
  
  return(result)
}

dataset.trim <- function (dataset, factor) {
  library(foreach)
  
  result <- list()
  
  for (i in names(dataset)) {
    if (is.list(dataset[[i]])) {
      result[[i]] <- foreach (j = dataset[[i]]) %do% {
        subset(j, subset = !is.na(factor))
      }
      names(result[[i]]) = names(dataset[[i]])
    } else {
      result[i] <- list(subset(dataset[i], subset = !is.na(factor)))
    }
  }
  
  return(result)
}

metric.individual <-
  function(df,
           grouping,
           value,
           factor,
           method = "wilcox",
           paired = FALSE,
           pal = "Dark2",
           lab.x = NULL,
           lab.legend = NULL) {
    library(vegan)
    library(ggplot2)
    library(ggpubr)
    library(foreach)
    
    result <- list()
    ##Perform statistical tests
    if (method == "t") {
      groups <- levels(factor(grouping))
      result$test <- foreach(i = groups) %do% {
        if (var.test(df[, which(colnames(df) == value)] ~ factor,
                     subset = grouping == i)$p.value > 0.05) {
          t.test(
            df[, which(colnames(df) == value)] ~ factor,
            subset = grouping == i,
            paired = paired,
            var.equal = T
          )
        } else {
          t.test(df[, which(colnames(df) == value)] ~ factor,
                 subset = grouping == i,
                 paired = paired)
        }
      }
      names(result$test) <- groups
    } else if (method == "wilcox") {
      groups <- levels(factor(grouping))
      result$test <- foreach(i = groups) %do% {
        wilcox.test(df[, which(colnames(df) == value)] ~ factor,
                    subset = grouping == i,
                    paired = paired)
      }
      names(result$test) <- groups
    }
    ##Prepare data for labels of significance
    test <-
      data.frame(
        grouping = groups,
        group1 = levels(factor)[1],
        group2 = levels(factor)[2],
        p = NA,
        y.position = max(df[, which(colnames(df) == value)]) * 1.02,
        p.signif = NA
      )
    for (i in 1:length(groups)) {
      test$p[i] = result$test[[i]]$p.value
      if (test$p[i] >= 0.05) {
        test$p.signif[i] = "ns"
      } else if (test$p[i] >= 0.01) {
        test$p.signif[i] = "*"
      } else if (test$p[i] >= 0.001) {
        test$p.signif[i] = "**"
      } else {
        test$p.signif[i] = "***"
      }
    }
    ##Plotting
    result["plot"] <- list(
      ggplot(df, aes(x = grouping, y = df[, which(colnames(df) == value)])) +
      scale_fill_brewer(palette = pal) +
      geom_boxplot(aes(fill = factor)) +
      labs(y = value, x = lab.x, fill = lab.legend) +
      stat_pvalue_manual(
        test,
        label = "p.signif",
        x = "grouping",
        hide.ns = T
      )
    )
    
    return(result)
  }

metric.overall <-
  function (df,
            value,
            factor,
            method = "wilcox",
            paired = FALSE,
            pal = "Dark2",
            lab.x = NULL,
            lab.legend = NULL) {
    library(vegan)
    library(ggplot2)
    library(ggpubr)
    
    result <- list()
    
    if (method == "t") {
      if (var.test(df[, which(colnames(df) == value)] ~ factor)$p.value > 0.05) {
        result["test"] <-
          list(t.test(df[, which(colnames(df) == value)] ~ factor,
                      var.equal = T,
                      paired = paired))
      } else {
        result["test"] <-
          list(t.test(df[, which(colnames(df) == value)] ~ factor,
                      paired = paired))
      }
    } else if (method == "wilcox") {
      result["test"] <-
        list(wilcox.test(df[, which(colnames(df) == value)] ~ factor,
                         paired = paired))
    } else if (method == "kruskal") {
      result["test"] <-
        list(kruskal.test(df[, which(colnames(df) == value)] ~ factor,
                          paired = paired))
    }
    test <-
      data.frame(
        group1 = levels(factor)[1],
        group2 = levels(factor)[2],
        p = result$test$p.value,
        y.position = max(df[, which(colnames(df) == value)]) * 1.02,
        p.signif = NA
      )
    if (test$p >= 0.05) {
      test$p.signif = "ns"
    } else if (test$p >= 0.01) {
      test$p.signif = "*"
    } else if (test$p >= 0.001) {
      test$p.signif = "**"
    } else {
      test$p.signif = "***"
    }
    
    result["plot"] <- list(
      ggplot(df,
             aes(x = factor,
                 y = df[, which(colnames(df) == value)])) +
      scale_fill_brewer(palette = pal) +
      geom_boxplot(aes(fill = factor)) +
      labs(y = value, x = lab.x, fill = lab.legend) +
      stat_pvalue_manual(test, label = "p.signif", hide.ns = T)
    )
    
    return(result)
  }

multi.adonis <- function(x, factors, sim.method, p.adjust.m)
{
  library(vegan)
  co = as.matrix(combn(unique(factors), 2))
  pairs = c()
  F.Model = c()
  R2 = c()
  p.value = c()
  
  for (elem in 1:ncol(co)) {
    ad = adonis(x[factors %in% c(as.character(co[1, elem]), as.character(co[2, elem])), ] ~ factors[factors %in% c(as.character(co[1, elem]), as.character(co[2, elem]))], method = sim.method)
    pairs = c(pairs, paste(co[1, elem], 'vs', co[2, elem]))
    F.Model = c(F.Model, ad$aov.tab[1, 4])
    R2 = c(R2, ad$aov.tab[1, 5])
    p.value = c(p.value, ad$aov.tab[1, 6])
  }
  p.adjusted <- p.adjust(p.value, method = p.adjust.m)
  pairw.res <- data.frame(pairs, F.Model, R2, p.value, p.adjusted)
  return(pairw.res)
}

multi.co.adonis <-
  function(x,
           covariant,
           factors,
           sim.method,
           p.adjust.m)
  {
    library(vegan)
    co = as.matrix(combn(unique(factors), 2))
    pairs = c()
    F.Model = c()
    R2 = c()
    p.value = c()
    
    for (elem in 1:ncol(co)) {
      ad = adonis(x[factors %in% c(as.character(co[1, elem]), as.character(co[2, elem])), ] ~ covariant[factors %in% c(as.character(co[1, elem]), as.character(co[2, elem]))] + factors[factors %in% c(as.character(co[1, elem]), as.character(co[2, elem]))], method = sim.method)
      pairs = c(pairs, paste(co[1, elem], 'vs', co[2, elem]))
      F.Model = c(F.Model, ad$aov.tab[2, 4])
      R2 = c(R2, ad$aov.tab[2, 5])
      p.value = c(p.value, ad$aov.tab[2, 6])
    }
    p.adjusted <- p.adjust(p.value, method = p.adjust.m)
    pairw.res <- data.frame(pairs, F.Model, R2, p.value, p.adjusted)
    return(pairw.res)
  }

multi.way.test <- function(df,
                           factor,
                           method = "kw",
                           trim = TRUE) {
  ##Prepare a blank data.frame for results
  result <-
    data.frame(
      "ID" = colnames(df),
      "p.value" = NA,
      "Adj.p.value" = NA,
      stringsAsFactors = FALSE
    )
  for (i in levels(factor)) {
    result[i] = NA
  }
  ##Perform the tests
  if (method == "kw") {
    for (i in 1:ncol(df)) {
      if (sum(df[, i]) == 0) {
        next()
      }
      temp <-
        kruskal.test(df[, i] ~ factor)
      result[i, 2] <- temp$p.value
      for (j in 4:ncol(result)) {
        result[i, j] = mean(subset(df[, i], subset = factor == levels(factor)[j - 3]))
      }
    }
  }
  
  if (trim) {
    result <- result[complete.cases(result[, 2]),]
  }
  result$Adj.p.value <- p.adjust(result$p.value, method = "fdr")
  
  return(result)
}

nmds.plot <-
  function(df,
           color,
           shape,
           method = "bray",
           shape_palette = c(10, 15),
           pal = "Dark2",
           conf = 0.95,
           alpha = 0.8,
           lab.shape = NULL,
           lab.color = NULL)
  {
    library(vegan)
    library(ggplot2)
    
    temp <- metaMDS(df)
    mds <- as.data.frame(temp$points)
    result <- ggplot(mds,
                     aes(
                       x = mds[, 1],
                       y = mds[, 2],
                       color = color,
                       shape = shape
                     )) +
      scale_shape_manual(values = shape_palette) +
      scale_color_brewer(palette = pal) +
      geom_point(size = 2, alpha = alpha) +
      stat_ellipse(level = conf) +
      labs(
        x = colnames(mds)[1],
        y = colnames(mds)[2],
        color = lab.color,
        shape = lab.shape,
        title = paste("Stress =", round(temp$stress, 4)),
        
      )
    return(result)
  }

pathway.cov <- function(ko_list) {
  library(foreach)
  library(doParallel)
  
  cl <- makeCluster(detectCores(), methods = F, useXDR = F)
  registerDoParallel(cl, cores = length(cl))
  result <-
    foreach(pathway = t(pathway_list), .combine = "c") %dopar% {
      ko <- unlist(strsplit(pathway[3], split = ","))
      sig_ko <- subset(ko, subset = ko %in% ko_list)
      cov <- length(sig_ko) / as.numeric(pathway[2])
    }
  stopCluster(cl)
  
  return(result)
}

pcoa.plot <-
  function(df,
           color,
           shape,
           method = "bray",
           pal = "Dark2",
           shape_palette = c(10, 15),
           conf = 0.95,
           alpha = 0.8,
           lab.shape = NULL,
           lab.color = NULL)
  {
    library(vegan)
    library(ggplot2)
    
    dm <- vegdist(df, method = method)
    pco_list <- cmdscale(dm, k = 2, eig = T)
    pco <- as.data.frame(pco_list$points)
    eig <- signif(pco_list$eig / sum(pco_list$eig) * 100, 2)
    colnames(pco) = paste("PCo", 1:2, ": ", eig[1:2], "%", sep = "")
    
    result <-
      ggplot(pco, aes(
        x = pco[, 1],
        y = pco[, 2],
        color = color,
        shape = shape
      )) +
      scale_shape_manual(values = shape_palette) +
      scale_color_brewer(palette = pal) +
      geom_point(size = 2, alpha = alpha) +
      stat_ellipse(level = conf) +
      labs(
        x = colnames(pco)[1],
        y = colnames(pco)[2],
        color = lab.color,
        shape = lab.shape
      )
    
    return(result)
  }

two.way.test <- function(df,
                         factor,
                         method = "wilcox",
                         alt = "two.sided",
                         paired = TRUE,
                         trim = TRUE) {
  
  ##Prepare a blank data.frame for results
  result <-
    data.frame(
      "ID" = colnames(df),
      "p.value" = NA,
      "Adj.p.value" = NA,
      "Fold.change" = NA,
      stringsAsFactors = FALSE
    )
  temp = paste(levels(factor), collapse = "/")
  colnames(result)[4] = paste0("Fold.change.", temp)
  ##Perform the tests
  ####ATTENTION: Permutation test (method = "perm") will ignore variables which has a prevalence less than 0.1 in either group.
  if (method == "perm") {
    library(foreach)
    library(doParallel)
    
    cl <- makeCluster(detectCores(), methods = F, useXDR = F)
    registerDoParallel(cl)
    test <- foreach (
      i = 1:ncol(df),
      .packages = c("foreach", "perm"),
      .errorhandling = "pass",
      .inorder = F
    ) %dopar% {
      temp <- foreach(j = levels(factor)) %do% {
        as.data.frame(proportions(table(subset(
          df[, i], factor == j
        ))))
      }
      foreach(k = temp) %do% {
        if (length(which(k$Var1 == 0)) != 0) {
          if (k[which(k$Var1 == 0),]$Freq > 0.9) {
            stop()
          }
        }
      }
      permTS(df[, i] ~ factor, alternative = alt)
    }
    stopCluster(cl)
    for (i in 1:length(test)) {
      if (is.null(test[[i]]$p.value)) {
        result[i, 2] <- NA
      } else {
        result[i, 2] <- test[[i]]$p.value
      }
      result[i, 4] <-
        log2(mean(subset(df[, i], factor == levels(factor)[1])) / mean(subset(df[, i], factor == levels(factor)[2])))
    }
  } else if (method == "t") {
    for (i in 1:ncol(df)) {
      if (sum(df[, i]) == 0) {
        next()
      }
      if (var.test(df[, i] ~ factor)$p.value > 0.05) {
        temp <- t.test(df[, i] ~ factor,
                       var.equal = T, paired = paired)
        result[i, 2] <- temp$p.value
        result[i, 4] <-
          log2(mean(subset(df[, i], factor == levels(factor)[1])) / mean(subset(df[, i], factor == levels(factor)[2])))
      } else {
        temp <- t.test(df[, i] ~ factor)
        result[i, 2] <- temp$p.value
        result[i, 4] <-
          log2(mean(subset(df[, i], factor == levels(factor)[1])) / mean(subset(df[, i], factor == levels(factor)[2])))
      }
    }
  } else if (method == "wilcox") {
    for (i in 1:ncol(df)) {
      if (sum(df[, i]) == 0) {
        next()
      }
      temp <-
        wilcox.test(df[, i] ~ factor, var.equal = T, paired = paired)
      result[i, 2] <- temp$p.value
      result[i, 4] <-
        log2(mean(subset(df[, i], factor == levels(factor)[1])) / mean(subset(df[, i], factor == levels(factor)[2])))
    }
  }
  
  if (trim) {
    result <- result[complete.cases(result[, 2]), ]
  }
  result$Adj.p.value <- p.adjust(result$p.value, method = "fdr")
  
  return(result)
}

unique.find <- function(df, factor, subset) {
  result <- data.frame("ID" = NA, "Mean" = 0, "Prevalence" = 0)
  
  for (i in 1:ncol(df)) {
    if (sum(subset(df[, i], subset = factor != subset)) == 0) {
      temp <- as.data.frame(proportions(table(subset(
        df[, i], factor == subset
      ))))
      result <- rbind(result, c(colnames(df)[i], mean(subset(df[, i], subset = factor == subset)), 1 - temp$Freq))
    }
  }
  result <- result[-1, ]
  
  return(result)
}

z.score.rep <- function(df, factor, cutoff = 0.95, cov = 0.5, pathway_list) {
  library(foreach)
  library(doParallel)
  
  temp_test <- two.way.test(df,
                            factor,
                            method = "perm",
                            alt = "greater")
  temp_test$Z.score <- qnorm(1 - temp_test$p.value)
  temp_test <- subset(temp_test, subset = temp_test$Z.score > -Inf)
  
  cl <- makeCluster(detectCores(), methods = F, useXDR = F)
  registerDoParallel(cl)
  temp_z <-
    foreach (i = 1:nrow(pathway_list), .combine = "c") %dopar% {
      sum(temp_test$Z.score[temp_test$ID %in% unlist(strsplit(pathway_list$ko[i], ","))]) / sqrt(pathway_list$ko_count[i])
    }
  names(temp_z) <- pathway_list$map
  temp_z.adj <-
    foreach(i = 1:nrow(pathway_list), .combine = "c") %dopar% {
      temp = 0
      set.seed(125)
      for (j in 1:1000) {
        temp[j] <-
          sum(sample(temp_test$Z.score, pathway_list$ko_count[i])) / sqrt(pathway_list$ko_count[i])
      }
      (temp_z[i] - mean(temp)) / sd(temp)
    }
  stopCluster(cl)
  
  temp_cov <- pathway.cov(colnames(df))
  result <- as.data.frame(cbind(temp_z, temp_z.adj, temp_cov))
  colnames(result) <- c("z.score", "z.score.adj", "pathway.cov")
  
  result <- subset(result, subset = pathway.cov > cov)
  result <- subset(result, subset = abs(result$z.score.adj) > qnorm(1 - (1 - cutoff) / 2))
  
  return(result)
}
