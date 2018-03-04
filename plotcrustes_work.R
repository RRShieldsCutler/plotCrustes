# Generate a more logical plot for procrustes comparisons
### By Robin Shields-Cutler
### Feburary 2018

library(ape)
library(vegan)
library(ggplot2)

# Read in beta diversity tables
unw_A <- read.delim('../unweighted_unifrac_asp_post_ampOTUnorm.txt',
                      header=1, row.names = 1, check.names = F)
unw_B <- read.delim('../unweighted_unifrac_stool_post_ampOTUnorm.txt',
                       header=1, row.names = 1, check.names = F)
# Metadata needs three columns - 
# 1. sample IDs (e.g. "P21_pre_stool")
# 2. the unifying participant/unit ID (e.g."P21") 
# 3. the binary metadata group, corresponding to the two distance input matrices (e.g. "pre" or "post")
meta <- read.delim('../plotcrust_metadata_aspstool.txt',
                   header=1, row.names=1, check.names = F)

thing <- 'patient'  # The meta header for the participant ID column
group <- 'group'  # The meta header for the metadata group category (which of the two distance matrices)
groups <- c('asp','stool')  # Names of the two categories in the group
metaA <- meta[meta[,group] == groups[1],]  # Split the metadata into the two groups
metaB <- meta[meta[,group] == groups[2],]

metaA <- metaA[order(metaA[,thing]),]  # Sort the metadata by unifying ID
metaB <- metaB[order(metaB[,thing]),]

# CRITICAL:
# Ensure that the original distance matrices are in the same order by participant
# Uses the order generated from the metadata dataframe
unw_A <- unw_A[rownames(metaA),rownames(metaA)]
unw_B <- unw_B[rownames(metaB),rownames(metaB)]

# Get the principal coordinates
pcoa_A <- pcoa(unw_A)$vectors
for(c in 1:ncol(pcoa_A)){
  colnames(pcoa_A)[c] <- paste0("PC",c)
}
pcoa_B <- pcoa(unw_B)$vectors
for(c in 1:ncol(pcoa_B)){
  colnames(pcoa_B)[c] <- paste0("PC",c)
}

crusty <- procrustes(pcoa_A, pcoa_B, symmetric = T)  # Run Procrustes
crust_test_p <- protest(pcoa_A, pcoa_B, permutations = how(nperm = 999))$signif
A_crust <- data.frame(crusty$X)  # Recover the first group's coordinates
B_crust <- data.frame(crusty$Yrot)  # Recover the second group's coordinates
colnames(B_crust) <- colnames(A_crust)
ncoords = as.numeric(ncol(A_crust))
A_crust <- merge(A_crust, metaA, by=0)
B_crust <- merge(B_crust, metaB, by=0)
rownames(A_crust) <- A_crust[,1]; A_crust[,1] <- NULL  # Merge makes rownames column
rownames(B_crust) <- B_crust[,1]; B_crust[,1] <- NULL
sample_ids <- A_crust[,thing]  # Get all the unifying participant IDs

real_dist <- data.frame(matrix(nrow = length(sample_ids), ncol = 3))  # Initialize the dataframe
colnames(real_dist) <- c('sampleID_A', 'sampleID_B', 'distance')
# Loop through each participant to get the multidimensional distance between their rotated points
for (i in 1:length(sample_ids)) {
  ix <- as.character(sample_ids[i])
  A_ix <- A_crust[A_crust[,thing] == ix, 1:ncoords]  # Keep all the PC axes
  B_ix <- B_crust[B_crust[,thing] == ix, 1:ncoords]
  AB_mat <- rbind(A_ix, B_ix)
  AB_dist <- matrix(dist(AB_mat, method = 'euclidean'))  # Calculates the distance
  real_dist[i,1] <- as.character(rownames(A_ix))  # Fill in the dataframe
  real_dist[i,2] <- as.character(rownames(B_ix))
  real_dist[i,3] <- as.numeric(AB_dist[1])
}

A_crust$timepoint <- 'A'; B_crust$timepoint <- 'B'
A_crust$realperm <- 'real'; B_crust$realperm <- 'real'

# The "Procrustes Distance"
# pro_dist <- sqrt(sum(real_dist$distance^2))


### Single permutation ###

metaBp <- metaB
# Permute the participant IDs to scramble the pair matching
# Later, people will be paired with non-self points for distance

metaBp[,thing] <- sample(x = metaBp[,thing], size = length(metaBp[,thing]), replace = F)


# Same as above
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

### Multiple permutations ###
### As used for the plots on the GitHub readme page

# Repeat the permutation "j" times
real_dist$real_perm <- 'true_distance'
dist_plot <- real_dist
PCOA_plot_many <- rbind(A_crust, B_crust)  # Start with the real data, then add the permutations
for (j in 1:4) { 
  metaBp <- metaB
  metaBp[,thing] <- sample(x = metaBp[,thing], size = length(metaBp[,thing]), replace = F)
  metaBp <- metaBp[order(metaBp[,thing]),]
  unw_Bp <- unw_B[rownames(metaBp),rownames(metaBp)]
  pcoa_Bp <- pcoa(unw_Bp)$vectors
  for(c in 1:ncol(pcoa_Bp)){
    colnames(pcoa_Bp)[c] <- paste0("PC",c)
  }
  crusty_p <- procrustes(pcoa_A, pcoa_Bp, symmetric = T)
  crust_test_perm_p <- protest(pcoa_A, pcoa_Bp, permutations = how(nperm = 999))$signif
  cat(paste0(crust_test_perm_p, '\n'))
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
  perm_dist_plot_2$real_perm <- paste0('permuted',j)  # Keep track of the permuted data
  dist_plot <- rbind(dist_plot, perm_dist_plot_2)
  
  # Add the permuted data to the existing runs
  Ap_crust$timepoint <- 'A'; Bp_crust$timepoint <- 'B'
  Ap_crust$realperm <- paste0('permuted',j); Bp_crust$realperm <- paste0('permuted',j)
  PCOA_plot_many <- rbind(PCOA_plot_many, Ap_crust, Bp_crust)
}


### Plot the results ###
# Boxplot with scatter
pbox <- ggplot(dist_plot, aes(x=real_perm, y=distance, group=real_perm)) +
  geom_boxplot(outlier.colour = 'white') + geom_jitter(width = 0.2) + theme_classic() +
  theme(axis.text.x = element_text(angle=45, hjust=1), axis.title.x = element_blank(),
        axis.text = element_text(color='black'))
# To save
ggsave(pbox, filename = '../distance_boxplots_aspstool_10x.png', width = 5, height = 3.5, dpi = 300)

# Set of typical procrustes plots
pscat <- ggplot(PCOA_plot_many) + geom_point(aes(x=PC1, y=PC2, color=timepoint)) +
  geom_line(aes(x=PC1, y=PC2, group=patient)) +
  theme_classic() + facet_grid(. ~ realperm) +
  theme(axis.text = element_text(color='black', size = 6))
# To save
ggsave(pscat, filename = '../permute10x_aspstool.png', width = 10, height = 1.5, dpi = 300)

