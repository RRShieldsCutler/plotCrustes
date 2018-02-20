# Generate a more logical plot for procrustes comparisons
### By Robin Shields-Cutler
### Feburary 2018

library(ape)
library(vegan)
library(ggplot2)

# Read in beta diversity tables
unw_A <- read.delim('../unweighted_unifrac_stool_pre_ampOTUnorm.txt',
                      header=1, row.names = 1, check.names = F)
unw_B <- read.delim('../unweighted_unifrac_stool_post_ampOTUnorm.txt',
                       header=1, row.names = 1, check.names = F)
meta <- read.delim('../plotcrust_metadata.txt',
                   header=1, row.names=1, check.names = F)

thing <- 'patient'
group <- 'group'
groups <- c('pre','post')
metaA <- meta[meta[,group] == groups[1],]
metaB <- meta[meta[,group] == groups[2],]

metaA <- metaA[order(metaA[,thing]),]
metaB <- metaB[order(metaB[,thing]),]

unw_A <- unw_A[rownames(metaA),rownames(metaA)]
unw_B <- unw_B[rownames(metaB),rownames(metaB)]

pcoa_A <- pcoa(unw_A)$vectors
for(c in 1:ncol(pcoa_A)){
  colnames(pcoa_A)[c] <- paste0("PC",c)
}

pcoa_B <- pcoa(unw_B)$vectors
for(c in 1:ncol(pcoa_B)){
  colnames(pcoa_B)[c] <- paste0("PC",c)
}

crusty <- procrustes(pcoa_A, pcoa_B, symmetric = T)
crust_test_p <- protest(pcoa_A, pcoa_B, permutations = how(nperm = 999))
A_crust <- data.frame(crusty$X)
B_crust <- data.frame(crusty$Yrot)
colnames(B_crust) <- colnames(A_crust)
ncoords = as.numeric(ncol(A_crust))
A_crust <- merge(A_crust, metaA, by=0)
B_crust <- merge(B_crust, metaB, by=0)
rownames(A_crust) <- A_crust[,1]; A_crust[,1] <- NULL
rownames(B_crust) <- B_crust[,1]; B_crust[,1] <- NULL
# B_crust$sample <- lapply(X = rownames(B_crust), FUN = function(xx)
#   strsplit(x=xx, split = 'endo', fixed = T)[[1]][1])
sample_ids <- A_crust[,thing]

real_dist <- data.frame(matrix(nrow = length(sample_ids), ncol = 3))
colnames(real_dist) <- c('sampleID_A', 'sampleID_B', 'distance')
for (i in 1:length(sample_ids)) {
  ix <- as.character(sample_ids[i])
  A_ix <- A_crust[A_crust[,thing] == ix, 1:ncoords]
  B_ix <- B_crust[B_crust[,thing] == ix, 1:ncoords]
  AB_mat <- rbind(A_ix, B_ix)
  AB_dist <- matrix(dist(AB_mat, method = 'euclidean'))
  real_dist[i,1] <- as.character(rownames(A_ix))
  real_dist[i,2] <- as.character(rownames(B_ix))
  real_dist[i,3] <- as.numeric(AB_dist[1])
}

A_crust$timepoint <- 'A'; B_crust$timepoint <- 'B'
A_crust$realperm <- 'real'; B_crust$realperm <- 'real'
pro_dist <- sqrt(sum(real_dist$distance^2))


#Permute

metaBp <- metaB
metaBp[,thing] <- sample(x = metaBp[,thing], size = length(metaBp[,thing]), replace = F)

crusty_p <- procrustes(pcoa_A, pcoa_B, symmetric = T)
crust_test_p <- protest(pcoa_A, pcoa_B, permutations = how(nperm = 99))$signif
Ap_crust <- data.frame(crusty_p$X)
Bp_crust <- data.frame(crusty_p$Yrot)
colnames(Bp_crust) <- colnames(Ap_crust)
Ap_crust <- merge(Ap_crust, metaA, by=0)
Bp_crust <- merge(Bp_crust, metaBp, by=0)
rownames(Ap_crust) <- Ap_crust[,1]; Ap_crust[,1] <- NULL
rownames(Bp_crust) <- Bp_crust[,1]; Bp_crust[,1] <- NULL
sample_ids <- Ap_crust[,thing]

perm_dist <- data.frame(matrix(nrow = length(sample_ids), ncol = 3))
colnames(perm_dist) <- c('sampleID_A', 'sampleID_B', 'distance')
for (i in 1:length(sample_ids)) {
  ix <- as.character(sample_ids[i])
  # ixp <- as.character(sample_id_perm[i])
  Ap_ix <- Ap_crust[Ap_crust[,thing] == ix, 1:14]
  Bp_ix <- Bp_crust[Bp_crust[,thing] == ix, 1:14]
  ABp_mat <- rbind(Ap_ix, Bp_ix)
  ABp_dist <- matrix(dist(ABp_mat, method = 'euclidean'))
  perm_dist[i,1] <- as.character(rownames(Ap_ix))
  perm_dist[i,2] <- as.character(rownames(Bp_ix))
  perm_dist[i,3] <- as.numeric(ABp_dist[1])
}

perm_pro_dist <- sqrt(sum(perm_dist$distance^2))

# perm_dist_plot <- perm_dist
# perm_dist_plot$real_perm <- 'permuted_distance'
# real_dist$real_perm <- 'true_distance'
# dist_plot <- rbind(real_dist, perm_dist_plot)

# ggplot(dist_plot, aes(x=real_perm, y=distance, fill=real_perm)) +
  # geom_boxplot() + theme_classic()


# Ap_crust$timepoint <- 'A'; Bp_crust$timepoint <- 'B'
# Ap_crust$realperm <- 'permuted'; Bp_crust$realperm <- 'permuted'
# 
# PCOA_plot <- rbind(A_crust, B_crust, Ap_crust, Bp_crust)
# PCOA_plot$SampleID <- rownames(PCOA_plot)
# 
# ggplot(PCOA_plot) + geom_point(aes(x=PC1, y=PC2, color=timepoint)) +
#   geom_line(aes(x=PC1, y=PC2, group=patient)) +
#   theme_classic() + facet_grid(. ~ realperm)

# make a bunch
real_dist$real_perm <- 'true_distance'
dist_plot <- real_dist
PCOA_plot_many <- rbind(A_crust, B_crust)
for (j in 1:9) { 
  metaBp <- metaB
  metaBp[,thing] <- sample(x = metaBp[,thing], size = length(metaBp[,thing]), replace = F)
  
  crusty_p <- procrustes(pcoa_A, pcoa_B, symmetric = T)
  crust_test_perm_p <- protest(pcoa_A, pcoa_B, permutations = how(nperm = 999))$signif
  Ap_crust <- data.frame(crusty_p$X)
  Bp_crust <- data.frame(crusty_p$Yrot)
  colnames(Bp_crust) <- colnames(Ap_crust)
  Ap_crust <- merge(Ap_crust, metaA, by=0)
  Bp_crust <- merge(Bp_crust, metaBp, by=0)
  rownames(Ap_crust) <- Ap_crust[,1]; Ap_crust[,1] <- NULL
  rownames(Bp_crust) <- Bp_crust[,1]; Bp_crust[,1] <- NULL
  sample_ids <- Ap_crust[,thing]
  
  perm_dist <- data.frame(matrix(nrow = length(sample_ids), ncol = 3))
  colnames(perm_dist) <- c('sampleID_A', 'sampleID_B', 'distance')
  for (i in 1:length(sample_ids)) {
    ix <- as.character(sample_ids[i])
    # ixp <- as.character(sample_id_perm[i])
    Ap_ix <- Ap_crust[Ap_crust[,thing] == ix, 1:ncoords]
    Bp_ix <- Bp_crust[Bp_crust[,thing] == ix, 1:ncoords]
    ABp_mat <- rbind(Ap_ix, Bp_ix)
    ABp_dist <- matrix(dist(ABp_mat, method = 'euclidean'))
    perm_dist[i,1] <- as.character(rownames(Ap_ix))
    perm_dist[i,2] <- as.character(rownames(Bp_ix))
    perm_dist[i,3] <- as.numeric(ABp_dist[1])
  }
  perm_dist_plot_2 <- perm_dist
  perm_dist_plot_2$real_perm <- paste0('permuted',j)
  dist_plot <- rbind(dist_plot, perm_dist_plot_2)
  
  Ap_crust$timepoint <- 'A'; Bp_crust$timepoint <- 'B'
  Ap_crust$realperm <- paste0('permuted',j); Bp_crust$realperm <- paste0('permuted',j)
  PCOA_plot_many <- rbind(PCOA_plot_many, Ap_crust, Bp_crust)
}

pbox <- ggplot(dist_plot, aes(x=real_perm, y=distance, group=real_perm)) +
  geom_boxplot(outlier.colour = 'white') + geom_jitter(width = 0.2) + theme_classic() +
  theme(axis.text.x = element_text(angle=45, hjust=1), axis.title.x = element_blank(),
        axis.text = element_text(color='black'))
ggsave(pbox, filename = '../distance_boxplots10x.png', width = 5, height = 3.5, dpi = 300)

pscat <- ggplot(PCOA_plot_many) + geom_point(aes(x=PC1, y=PC2, color=timepoint)) +
  geom_line(aes(x=PC1, y=PC2, group=patient)) +
  theme_classic() + facet_grid(. ~ realperm) +
  theme(axis.text = element_text(color='black', size = 6))
ggsave(pscat, filename = '../permute10x.png', width = 10, height = 1.5, dpi = 300)

