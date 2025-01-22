# SET UP ----

setwd("~/Library/CloudStorage/OneDrive-udd.cl/20240823_final_discordance/02_discordance")
source("~/Library/CloudStorage/GoogleDrive-berebolledo@gmail.com/My Drive/UDD/research/bioinformatics/SABIO/projects/03_github_repos/fondecyt_11220642/05_discordance-functions.R")

library(readxl)
library(dplyr)
library(stringr)
library(readr)
library(ggpubr)
library(gridExtra)

set.seed(11220642)

# MITOCARTA ----    

mitopath <- read_excel("inputs/Human.MitoCarta3.0.xls", 
                       sheet = "C MitoPathways", 
                       col_types = c("skip", "text", "text", "text"))

mitopath <- mitopath[complete.cases(mitopath),]

# Generate gene lists

oxphos <- getGenes(mitopath, "OXPHOS")

high <- read_delim("inputs/genes.high-mt.bed", 
                   delim = "\t", escape_double = FALSE, 
                   col_names = FALSE, trim_ws = TRUE)
high <- high$X4
low <- read_delim("inputs/genes.low-mt.bed", 
                  delim = "\t", escape_double = FALSE, 
                  col_names = FALSE, trim_ws = TRUE)
low <- low$X4
allmitogenes <- c()
for (path in mitopath$MitoPathway){
  allmitogenes <- c(allmitogenes, getGenes(mitopath, path))
}
allmitogenes <- unique(allmitogenes)

# Change accordingly, possible values:

# mygenelist <- oxphos
# q22deltafile <- "22q/22q_delta_oxphos.txt" 

# mygenelist <- high
# q22deltafile <- "22q/22q_delta_high.txt" 

# mygenelist <- low
# q22deltafile <- "22q/22q_delta_low.txt" 

mygenelist <- allmitogenes
q22deltafile <- "22q/22q_delta_allmito.txt"


# 1000 genomes megatable and delta MND ----
an1 <- read.delim("g1k/g1k.an1.discordance.tsv")
nonmitogenes <- an1$gene[!an1$gene %in% allmitogenes] 
population <- read.csv("g1k/20130606_g1k_3202_samples_ped_population.txt", sep="")
population <- population[,c("SampleID", "Sex", "Population", "Superpopulation")]
haplogroups <- read.csv("g1k/mito_anc_3202_short_version.txt", sep = "\t", quote = "")
haplogroups$macroHap <- substr(haplogroups$Haplogroup,1,1)
haplogroups <- haplogroups[,c("SampleID", "macroHap")]
pedigree <- read.csv("g1k/1kGP.3202_samples.pedigree_info.txt", sep="")
pedigree <- pedigree[,1:3]
colnames(pedigree) <- c("SampleID", "FatherID", "MotherID")

check1 <- checkNA(an1)
if (length(check1$narows) > 0 ){
  an1 <- an1[-check1$narows,]
}

allmito1 <- create_tables(an1, mygenelist)
g1kdelta <- allmito1$deltaMother[complete.cases(allmito1$deltaMother)]
g1kdelta <- data.frame(delta = g1kdelta, label = rep("1kGP-pairs", length(g1kdelta)))

# Figure 1A: Population specific MND distribution ----

g1k <- read.delim("g1k/g1k.an1.discordance.tsv")
sg <- read.delim("sg/sg.an1.discordance.tsv")
sgfam <- sg[,colnames(sg)[!grepl("F", colnames(sg))]]
sg <- sg[,colnames(sg)[!grepl("P", colnames(sg))]]

merged <- merge(g1k,sg)
merged <- merged[merged$gene %in% mygenelist,] #CHANGE GENE LIST <---
colMeans <- apply(merged, 2, FUN = function(x) mean(as.numeric(x),na.rm = T))

df <- data.frame(SampleID = colnames(merged), meanMND = colMeans)
df <- df[2:nrow(df),]

population <- read.delim("g1k/20130606_g1k_3202_samples_ped_population.txt", sep = " ")

pop2 = data.frame(FamilyID = colnames(sg)[2:ncol(sg)],
                  SampleID = colnames(sg)[2:ncol(sg)],
                  FatherID = rep(0, 63),
                  MotherID = rep(0, 63),
                  Sex      = rep(0, 63),
                  Population = rep("CHL", 63),
                  Superpopulation = rep("AMR", 63)
)

population <- rbind(population, pop2)
supertable <- merge(df, population)


ggboxplot(supertable, x = "Population", y = "meanMND",
          color = "Population", lwd = 1, fatten = 0.5,
          add.params = list(color = "Population", size = 1.5, shape = 19),
          add = "jitter") +
  labs(title = "", x= "Population", y = "Mean MND") +
  theme(legend.position = "none") +
  stat_compare_means(method="anova", label.y = 1.2) 

# Black and white

a <- ggboxplot(supertable, x = "Population", y = "meanMND",
          lwd = 1, fatten = 0.5,
          add.params = list(alpha = 0.5, size = 1.5, shape = 19),
          add = "jitter") +
  labs(title = "", x= "Population", y = "Mean MND") +
  theme(legend.position = "none") +
  stat_compare_means(method="anova", label.y = 1.5) +
  scale_y_continuous("Mean MND", breaks = c(0, 0.5, 1), labels = c(0, 0.5, 1), limits = c(0, 2))
a
summary(aov(meanMND ~ Population, data = supertable))
summary(aov(meanMND ~ Population, data = supertable))[[1]][["Pr(>F)"]][1]


# Correlation with european ancestry
eur_g1k <- read.table("../04_global_ancestry/trios.eur_ancestry.tsv", sep = "\t", header = TRUE)
colnames(eur_g1k) <- c("SampleID", "EUR_ANCESTRY")
eur_sg <- read.table("../04_global_ancestry/sg.eur_ancestry.tsv", sep = "\t", header = TRUE)
colnames(eur_sg) <- c("SampleID", "EUR_ANCESTRY")
eur_sg$SampleID <- gsub("-", ".", eur_sg$SampleID)
eur_ancestry <- rbind(eur_g1k, eur_sg)
supertable <- merge(supertable,eur_ancestry)
ggscatter(supertable, x = "meanMND", y = "EUR_ANCESTRY",
          add = "reg.line", shape = 19, alpha = 0.5,
          add.params = list(color = "blue", fill = "lightgray"), 
          conf.int = TRUE) +
  stat_cor(method = "pearson", label.x = 0, label.y = 0.5, r.accuracy = 0.001)+
  labs(x = "mean MND", y = "Global European ancestry") +
  scale_y_continuous("Global European ancestry", breaks = c(0, 0.25, 0.5, 0.75, 1), labels = c(0.00, 0.25, 0.50, 0.75, 1.00), limits = c(0, 1))+
  scale_x_continuous("mean MND", breaks = c(0, 0.25, 0.5, 0.75, 1), labels = c(0.00, 0.25, 0.50, 0.75, 1.00), limits = c(0, 1))

cor.test(supertable$meanMND, supertable$EUR_ANCESTRY)$p.value

supertable %>%
  group_by(Population) %>%
  summarise(avg = mean(meanMND), sd = sd(meanMND))


# Figure 1B: Haplogroup specific MND distribution ----

hapg1k <- read.delim("g1k/mito_anc_3202_short_version.txt")
hapg1k <- hapg1k[, 1:2]
hapsg <- read.delim("sg/sg.haplogroups.tsv")
hapsg <- hapsg[, 1:2]
haps <- rbind(hapg1k, hapsg)

haps$macroHap <- substr(haps$Haplogroup,1,1)
supertable$SampleID <- gsub("\\.", "-", supertable$SampleID)
supertable2 <- merge(supertable, haps)
supertable2 <- supertable2[supertable2$macroHap %in% c("A", "B", "C", "D"),]

b <- plothap3(supertable2)
b
summary(aov(meanMND ~ macroHap, data = supertable2))


supertable2 %>%
  group_by(macroHap) %>%
  summarise(avg = mean(meanMND), sd = sd(meanMND), var = var(meanMND))

ggarrange(a,b, ncol = 2, labels = "AUTO")
# 700 x 385

# Figure 2: Cohort specific MND ----

sg <- read.delim("sg/sg.an1.discordance.tsv")
sg <- sg[,colnames(sg)[!grepl("P", colnames(sg))]]

so <- read.delim("sophia/sophia.an1.discordance.tsv")
nw <- read.delim("nwu/nwu.an1.discordance.tsv")
q2 <- read.delim("22q/22q.an1.discordance.tsv")

merge1 <- merge(sg,so)
merge2 <- merge(merge1, nw)
merge3 <- merge(merge2, q2)
merged <- merge3
merged <- merged[merged$gene %in% mygenelist,] #CHANGE GENE LIST <---

colMeans <- apply(merged, 2, FUN = function(x) mean(as.numeric(x),na.rm = T))

df <- data.frame(SampleID = colnames(merged), meanMND = colMeans)
df <- df[2:nrow(df),]

labels <- c(rep("DRD-healthy",   length(colnames(sg)) - 1),
            rep("DRD-affected", length(colnames(so)) - 1),
            rep("22q-CHL", length(colnames(nw)) - 1),
            rep("22q-ARG", length(colnames(q2)) - 1)
)

df$Population <- labels


my_comparisons <- list(c("DRD-healthy", "DRD-affected"), c("DRD-healthy", "22q-CHL"), c("DRD-healthy", "22q-ARG"),
                       c("DRD-affected", "22q-CHL"), c("DRD-affected", "22q-ARG"),
                       c("22q-CHL","22q-ARG"))


d <- ggboxplot(df, x = "Population", y = "meanMND",
          #color = "Population", 
          lwd = 1, fatten = 0.5,
          #add.params = list(color = "Population", size = 1.5, shape = 19),
          add.params = list( size = 1.5, shape = 19, alpha = 0.5),
          add = "jitter") +
  labs(title = "", x= "Cohort", y = "mean MND") +
  theme(legend.position = "none") +
  stat_compare_means(method = "t.test", comparisons = my_comparisons, label = "p.signif")  +
  stat_compare_means(method="anova", label.y = 1.6, label.x = 0.75) +
  scale_y_continuous("mean MND", breaks = c(0, 0.5, 1), labels = c(0.0, 0.50, 1.00), limits = c(0, 1.7))
d
summary(aov(meanMND ~ Population, data = df))

#allmitodf <- df
#oxphosdf <- df
#highdf <- df
#lowdf <- df

#allmitodf$class = "mitogenes"
#oxphosdf$class = "oxphos"
#highdf$class = "high-mt"
#lowdf$class = "low-mt"

#bigdf <- rbind(allmitodf, oxphosdf, highdf, lowdf)

#bigdf$class <- factor(bigdf$class, levels = c("mitogenes", "oxphos", "high-mt", "low-mt"))

#ggboxplot(bigdf, x = "Population", y = "meanMND", facet.by = "class",
#          lwd = 1, fatten = 0.5,
#          add.params = list( size = 1.5, shape = 19, alpha = 0.5),
#          add = "jitter") +
#  labs(title = "", x= "", y = "mean MND") +
#  theme(legend.position = "none") +
#  stat_compare_means(method = "t.test", comparisons = my_comparisons, label = "p.signif")  +
#  stat_compare_means(method="anova", label.y = 1.6, label.x = 0.75) +
#  scale_y_continuous("mean MND", breaks = c(0, 0.5, 1), labels = c(0.0, 0.50, 1.00), limits = c(0, 1.7))
  



df %>%
  group_by(Population) %>%
  summarise(mean = mean(meanMND), sd = sd(meanMND), var =var(meanMND))

aggregate(df$meanMND, by = list(df$Population), FUN = function(x) shapiro.test(x)$p.value)


var.test(df[df$Population == "22q-ARG",]$meanMND,df[df$Population == "22q-CHL",]$meanMND)
var(df[df$Population == "22q-ARG",]$meanMND) / var(df[df$Population == "22q-CHL",]$meanMND)

#ggarrange(c,d, ncol = 2, labels = "AUTO")
# 950 x 490

# Figure 3: Delta comparisons ----

sg <- read.delim("sg/sg.an1.discordance.tsv")
sgfam <- sg[,colnames(sg)[!grepl("F", colnames(sg))]]

sgfam <-  sgfam[sgfam$gene %in% mygenelist,] # CHANGE GENE LIST <---
colMeans <- apply(sgfam, 2, FUN = function(x) mean(as.numeric(x),na.rm = T))
sgfam <- data.frame(SampleID = colnames(sgfam), meanMND = colMeans)
sgfam <- sgfam[2:nrow(sgfam),]
probands <- sgfam[grepl("P", sgfam$SampleID),]
probands <- probands[!grepl("157", probands$SampleID),]
mothers  <- sgfam[grepl("M", sgfam$SampleID),]

sgdelta <- cbind(probands, mothers)
sgdelta$delta <- (sgdelta[,2] - sgdelta[,4])/sgdelta[,4]
sgdelta$label <- rep("DRD-pairs", nrow(sgdelta))


q2 <- read.delim("22q/22q.an1.discordance.tsv")
q2 <- q2[q2$gene %in% mygenelist,] # CHANGE GENE LIST <---
colMeans <- apply(q2, 2, FUN = function(x) mean(as.numeric(x),na.rm = T))
q2 <- data.frame(SampleID = colnames(q2), meanMND = colMeans)
q2 <- q2[2:nrow(q2),]

q2ped <- read.delim("22q/22q_pedigree.txt")
q2ped <- q2ped[q2ped$MotherID !=0,]

q2delta <- read.delim(q22deltafile, header = FALSE) #CHANGE GENE LIST <---
q2delta$delta <- (q2delta$V2 - q2delta$V4)/q2delta$V4
q2delta$label <- rep("22q-ARG-pairs", nrow(q2delta))

comparedelta <- rbind(g1kdelta, sgdelta[,c(5,6)], q2delta[, c(5,6)])

ggboxplot(comparedelta, x = "label", y = "delta",
          #color = "label", 
          lwd = 1, fatten = 0.5,
          #add.params = list(color = "label", size = 1.5, shape = 19),
          add.params = list(size = 1.5, shape = 19, alpha = 0.5),
          add = "jitter") +
  labs(title = "", x= "Cohort", y = "∆MND") +
  theme(legend.position = "none") +
  stat_compare_means(method="kruskal", label.y = 1, label.x = 2.5) +
  scale_y_continuous("∆MND", breaks = c(-1, -0.5, 0, 0.5, 1), labels = c(-1, -0.5, 0, 0.5, 1), limits = c(-1, 1))
  

comparedelta %>%
  group_by(label) %>%
  summarise(median = median(delta), min = min(delta), max = max(delta), var = var(delta))

var.test(comparedelta[comparedelta$label == "1kGP-pairs",]$delta,comparedelta[comparedelta$label == "DRD-pairs",]$delta)
var.test(comparedelta[comparedelta$label == "1kGP-pairs",]$delta,comparedelta[comparedelta$label == "22q-ARG-pairs",]$delta)


x <- seq(-1, 1, length = 80)
hist(comparedelta$delta, col = "white", border = "white", prob = TRUE, 
     ylim = c(0, 5), xlim=c(-1, 1), yaxt = "n", ylab = "", xlab = "∆MND",
     main = "")

fun1 <- dnorm(x, mean = mean(g1kdelta$delta), sd = sd(g1kdelta$delta))
lines(x, fun1, col = "tomato" , lwd = 4)

fun2 <- dnorm(x, mean = mean(sgdelta$delta), sd = sd(sgdelta$delta))
lines(x, fun2, col = "forestgreen", lwd = 4)

fun3 <- dnorm(x, mean = mean(q2delta$delta), sd = sd(q2delta$delta))
lines(x, fun3, col = "gold", lwd = 4)




# 500 x 364
