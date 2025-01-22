# Calculate MND and other functions ----

checkNA <- function(dataset) {
  
  out <- list()
  
  out$rowsummary <- summary(apply(dataset[,2:ncol(dataset)], 2, function(x) sum(is.na(x))))
  
  narows <- c()
  for (col in colnames(dataset)){
    narows <- c(narows, which(is.na(dataset[[col]])))
  }
  
  out$colsummary <- summary(data.frame(table(narows))$Freq)
  out$uniquerows <- length(unique(narows))
  
  out$narows <- unique(narows)
  
  return(out)
  
}

getGenes<- function(mitocarta, pathway){
  
  genes.tmp <- (mitocarta %>% 
                  filter(MitoPathway == pathway) %>% 
                  select(Genes))$Genes
  
  gene.list <- strsplit(genes.tmp, split = ",")
  gene.list <- unlist(lapply(gene.list, FUN = str_trim))
  return(gene.list)
  
}

genomeDiff <- function(data, genelist, sample, alt = "two.sided"){
  
  
  
  expectedDiscordance <-c()
  for (i in 1:1000){
    
    sampledRows <- sample(1:nrow(data), length(genelist))
    meanRandomDiscordance <- mean(data[sampledRows, sample])
    expectedDiscordance <- c(expectedDiscordance, meanRandomDiscordance)
  }
  
  mitoDiscordance <- data[data$gene %in% genelist, sample]
  
  obsDiscordance <- mean(mitoDiscordance)
  expDiscordance <- mean(expectedDiscordance)
  
  testRes <- try(t.test(x = mitoDiscordance, mu = expDiscordance, alternative = alt), silent = TRUE)
  
  if (is(testRes, "try-error")) {
    return(list(sample = sample, ngenes = length(genelist), obsMND = obsDiscordance, expMND = expDiscordance, p = NA))
  } else {
    return(list(sample = sample, ngenes = length(genelist), obsMND = obsDiscordance, expMND = expDiscordance, p = testRes$p.value))
  }
  
}

calcDiff <- function(data, genelist){
  
  p.vals <- c()    
  for (sample in colnames(data)[2:ncol(data)]){
    p.vals <- c(p.vals, genomeDiff(data, genelist, sample)$p)
  }
  
  p.vals.adjusted <- p.adjust(p.vals, method = "fdr")
  results <- colnames(data)[which(p.vals.adjusted < 0.05) + 1]
  return(p.vals.adjusted )
}

calcObsMND <- function(data, genelist){
  
  obs.MND <- c()    
  for (sample in colnames(data)[2:ncol(data)]){
    obs.MND <- c(obs.MND, genomeDiff(data, genelist, sample)$obsMND)
  }
  
  return(obs.MND)
}

calcExpMND <- function(data, genelist){
  
  exp.MND <- c()    
  for (sample in colnames(data)[2:ncol(data)]){
    exp.MND <- c(exp.MND, genomeDiff(data, genelist, sample)$expMND)
  }
  
  return(exp.MND)
}        

bigTable <- function(ancestry, geneList){
  
  res <- data.frame(SampleID = colnames(ancestry)[2:ncol(ancestry)], 
                    obsMND = calcObsMND(ancestry, geneList), 
                    expMND = calcExpMND(ancestry, geneList), 
                    fdr = calcDiff(ancestry, geneList))
  df1 <- merge(res, population, by = "SampleID")
  df2 <- merge(df1, haplogroups, by = "SampleID")
  df3 <- merge(df2, pedigree, by = "SampleID" )
  
  
  delta <- c()
  deltafa <- c()
  motherObs <- c()
  fatherObs <- c()
  
  for (i in 1:nrow(df3)) {
    
    mother <- df3[i,]$MotherID
    
    if( mother != "0") {
      
      motherObsMND <- df3[which(df3$SampleID == mother), "obsMND"]
      
      if (length(motherObsMND) != 0 && motherObsMND != 0 ) {
        d <- (df3[i,]$obsMND - motherObsMND)/motherObsMND
        motherObsMND <- motherObsMND
        
      } else {
        d <- NA
        motherObsMND <- NA}
      
    } else {
      d <- NA
      motherObsMND <- NA
      
    }
    
    delta <- c(delta, d)
    motherObs <- c(motherObs, motherObsMND )
    
  }
  
  
  for (i in 1:nrow(df3)) {
    
    father <- df3[i,]$FatherID
    
    if( father != "0") {
      
      fatherObsMND <- df3[which(df3$SampleID == father), "obsMND"]
      
      if (length(fatherObsMND) != 0 && fatherObsMND != 0 ) {
        dfa <- (df3[i,]$obsMND - fatherObsMND)/fatherObsMND
        fatherObsMND <- fatherObsMND
        
      } else {
        dfa <- NA
        fatherObsMND <- NA}
      
    } else {
      dfa <- NA
      fatherObsMND <- NA
      
    }
    
    deltafa <- c(deltafa, dfa)
    fatherObs <- c(fatherObs, fatherObsMND )
    
  }
  
  df3$deltaMother <- delta
  df3$motherObsMND <- motherObs
  df3$deltaFather <- deltafa
  df3$fatherObsMND <- fatherObs
  
  
  return(df3)
}

create_tables <- function(ancestry,genelist){
  df <- bigTable(ancestry, genelist)
  df <- df[df$macroHap %in% c("A", "B", "C", "D"),]
  df$Population <- as.factor(df$Population)
  df$prop <- df$obsMND/df$expMND
  df$macroHap <- as.factor(df$macroHap)
  df$Sex <- factor(df$Sex, levels = c(1,2), labels = c("Men", "Women"))
  return(df)
}

# This function calculates a permutation test on the Mann-Whitney U test statistic

boot_test_2 <- function(concordant, discordant){
  
  set.seed(2023)
  
  ncon <- length(concordant)
  ndis <- length(discordant)
  
  obs <- wilcox.test(concordant, discordant)$statistic
  
  combined <- c(concordant, discordant)
  
  boot = 0
  for (i in 1:1000){
    
    newsample <- sample(combined, ncon+ndis, replace=FALSE)
    newcon <- newsample[1:ncon]
    newdis <- newsample[(ncon+1):(ncon+ndis)]
    newstat <- wilcox.test(newcon, newdis)$statistic
    #newstat <- t.test(newcon, newdis)$statistic
    #boot <- c(boot, newstat)
    if (newstat>=obs){
      boot <- boot + 1
    } else {
      boot <- boot + 0
    }
    
  }
  
  return(boot/1000) 
}

plotpop <- function(df){
  
  my_comparisons <- list(c("CLM", "MXL"),  c("CLM", "PEL"),  c("CLM", "PUR"),
                         c("MXL", "PEL"),  c("MXL", "PUR"),
                         c("PEL", "PUR"))
  
  ggboxplot(df, x = "Population", y = "obsMND",
            color = "Population", lwd = 1, fatten = 0.5,
            add.params = list(color = "Population", size = 1.5, shape = 19),
            add = "jitter") +
    labs(title = "", x= "Population", y = "Observed MND") +
    theme(legend.position = "none") +
    stat_compare_means(comparisons = my_comparisons, label = "p.signif")  
}

plothap <- function(df){
  
  my_comparisons <- list(c("A","B"), c("A","C"), c("A", "D"),
                         c("B", "C"), c("B", "D"),
                         c("C", "D"))
  
  ggboxplot(df, x = "macroHap", y = "obsMND",
            color = "gray60", lwd = 1, fatten = 0.5,
            add.params = list(color = "Population", size = 1.5, shape = 19),
            add = "jitter") +
    labs(title = "", x= "Macrohaplogroup", y = "Observed MND") +
    theme(legend.position = "none") +
    stat_compare_means(comparisons = my_comparisons, label = "p.signif")  
  
  
}

plotinter <- function(df1, df2){
  d1 <- df1[, c("Sex", "deltaMother")]
  d1 <- d1[complete.cases(d1),]
  d1$Haplotype <- rep("Haplotype 1", nrow(d1))
  
  d2 <- df2[, c("Sex", "deltaMother")]
  d2 <- d2[complete.cases(d2),]
  d2$Haplotype <- rep("Haplotype 2", nrow(d2))
  
  dat <- rbind(d1,d2)
  dat$Sex <- as.factor(dat$Sex)
  dat$Haplotype <- as.factor(dat$Haplotype)
  
  p <- ggboxplot(dat, x = "Sex", y = "deltaMother", facet.by = "Haplotype",
                 color = "Sex", lwd = 0.5, fatten = 0.5,
                 add.params = list(color = "Sex", size = 0.5, shape = 19),
                 add = "jitter") +
    labs(title = "", x= "", y = "∆MND") +
    theme(legend.position = "none")
  
  my_comparisons <- list( c("Men", "Women"))
  
  p <- p + stat_compare_means(comparisons = my_comparisons, method = "wilcox.test",
                              label = "p.signif", bracket.size = 0.7, size = 5,
                              label.y = 1)
  return(p)
}

plothap2 <- function(df){
  
  my_comparisons <- list(c("A","B"), c("A","C"), c("A", "D"),
                         c("B", "C"), c("B", "D"),
                         c("C", "D"))
  
  ggboxplot(df, x = "macroHap", y = "meanMND",
            color = "gray60", lwd = 1, fatten = 0.5,
            add.params = list(color = "Population", size = 1.5, shape = 19),
            add = "jitter") +
    labs(title = "", x= "Macrohaplogroup", y = "Mean MND") +
    theme(legend.position = "none") +
    stat_compare_means(comparisons = my_comparisons, label = "p.signif") +
    stat_compare_means(method = "kruskal", label.y = 1.5)
  
  
}

plothap3 <- function(df){
  
  my_comparisons <- list(c("A","B"), c("A","C"), c("A", "D"),
                         c("B", "C"), c("B", "D"),
                         c("C", "D"))
  
  ggboxplot(df, x = "macroHap", y = "meanMND",
            lwd = 1, fatten = 0.5,
            add.params = list(alpha = 0.5, size = 1.5, shape = 19),
            add = "jitter") +
    labs(title = "", x= "Macrohaplogroup", y = "Mean MND") +
    theme(legend.position = "none") +
    stat_compare_means(comparisons = my_comparisons, label = "p.signif") +
    stat_compare_means(method = "anova", label.y = 1.5) +
    scale_y_continuous("mean MND", breaks = c(0,0.5,1), labels = c(0,0.5,1), limits = c(0, 2))
  
  
  
}

