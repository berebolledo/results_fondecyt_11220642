library(ggplot2)

read_anc <- function(cohort, ancestry, acronym){
  df <- read.table(paste0(cohort,".",ancestry,"_ancestry.tsv"), header = T)
  colnames(df) <- c("sample", "estimate")
  df$anc <- rep(acronym, nrow(df))
  return(df)
}

dataset <- "1kGP-AMR" # 22q-ARG; 22q-CHL; DRD-healthy; DRD-affected

df <- rbind(read_anc(dataset, "afr", "AFR"),
      read_anc(dataset, "amr", "NAT"),
      read_anc(dataset, "eur", "EUR"))

df2 <- df[df$anc == "NAT",]
rownames(df2) <- df2$sample
df$sample <- factor(df$sample, levels = rownames(df2)[order(df2$estimate)] )

# Stacked
ggplot(df, aes(fill=anc, y=estimate, x=sample)) + 
  geom_bar(position="stack", stat="identity", width=1) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +
  ggtitle(dataset) + ylab("Global ancestry estimate")





