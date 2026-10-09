## Estimation of assemblies divergence
## R script

library(data.table)

files <- list.files(
  pattern = "^E_tef_pangenome\\..*\\.mash_triangle\\.txt$",
  full.names = TRUE
)

read_triangle <- function(file) {
  lines <- readLines(file)
  #lines <- readLines(files)

  n <- length(lines)
  
  rownms <- sapply(strsplit(lines, "\\s+"), `[`, 1)
  
  L <- matrix(0, n, n)
  rownames(L) <- rownms
  colnames(L) <- rownms
  
  for(i in seq_len(n)) {#i=2
    parts <- strsplit(lines[i], "\\s+")[[1]]
    l_parts<-length(parts)
    if (l_parts > 1){
      vals <- as.numeric(parts[-1])  # remove row name
      L[i, 1:(i-1)] <- vals
    } 
  }
  
  L
}

matrices <- lapply(files, read_triangle)

# Name list elements according to file names
names(matrices) <- basename(files)
array<-simplify2array(matrices)
dim(array)


mean_matrix <- apply(array,c(1, 2),mean)
dim(mean_matrix)
write.table(mean_matrix, "mean_divergence.txt", quote=F, sep="\t", 
            row.names=T, col.names=F)
#check
v<-rep(0,dim(array)[3])
for (i in 1 : dim(array)[3]){#i=5
  v[i]<-array[16,10,i]
  array[,,i][upper.tri(array[,,i])] <- t(array[,,i])[upper.tri(array[,,i])]
} 
mean(v)
mean_matrix[16,10]
max(mean_matrix)


###############
#plot per chromosome similarities
##############
#which chromosome shows bigger variation?

library(ggplot2)
library(tidyr)
library(viridis)
#need a df with a col per chromosome and rows that are the pairwise comparisons
coln<- c(colnames(matrices[[1]]))
#coln_new<-sub("_[^_]+kb_", "_", coln)
coln_new<-sub("_.*", "", coln)
rown<- c(rownames(matrices[[1]]))
rown_new<-sub("_.*", "", rown)
#df_rown<-apply(combn(rown_new, 2), 2, paste, collapse = "_")
df_rown<-apply(expand.grid(rown_new, coln_new), 1, paste0, collapse = "_")


df<- matrix(data=NA,nrow=length(df_rown),ncol=length(matrices))
rownames(df)<-df_rown
colnames(df)<-c('1A','1B','2A','2B','3A','3B','4A','4B','5A','5B','6A','6B','7A','7B','8A','8B','9A','9B','10A','10B')
for (i in 1:length(matrices)){#i=1
  matrices[[i]][upper.tri(matrices[[i]])] <- t(matrices[[i]])[upper.tri(matrices[[i]])]
  df[,i]<-c(matrices[[i]])
}
df<-as.data.frame(df)
which(df[,1] == max(df[,1]))

max<-c()
for (i in 1:ncol(df)){
  max[i]<-max(df[,i])
}

write.table(df, "divergence_bychr.txt", quote=F, sep="\t", 
            row.names=T, col.names=F)




df_long <- pivot_longer(
  df,
  cols = everything(),
  names_to = "Variable",
  values_to = "Value"
)

df_long$chr<-grep('A',df_long)
df_long$chr<-df_long$Variable
df_long$chr<-gsub('A','',df_long$chr)
df_long$chr<-gsub('B','',df_long$chr)
levels_order <- as.vector(rbind(paste0(1:10, "A"), paste0(1:10, "B")))
#levels_order <- c(1:10)
df_long$chr<-as.factor(df_long$chr)
df_long$Variable<-as.factor(df_long$Variable)
class(df_long$chr)
df_long<-as.data.frame(df_long)

#df_long$haplo_chr<-as.factor(c(1,1,2,2,3,3,4,4,5,5,6,6,7,7,8,8,9,9,10,10))
df_long$Variable<- factor(df_long$Variable, levels = levels_order)
levels_chr <- c(1:10)
df_long$chr<- factor(df_long$chr, levels = levels_chr)

class(df_long$Variable)



pdf('div_bychr.pdf', width=13, height=7)
ggplot(df_long, aes(x = Variable, y = Value, fill= chr)) +
  geom_boxplot() +
  #scale_fill_brewer(palette="BuPu")+
  #scale_fill_brewer(palette='Set1')+
  scale_fill_manual(values = c('#c58b3e','#72c977','#883880','#c01933','#6582d8','#bc495f','#d397ff','#43d2c2','#bc5b37','#639337'))+ 
  theme_bw() +
  labs( x = "Chromosome", y = "Divergence") +
  labs(fill = "Homeoulogous\nChromosome") +
  theme(
  # Axis titles
  axis.title.x = element_text(size = 20), #face = "bold"),
  axis.title.y = element_text(size = 20), #face = "bold"),
  #axis tick labels
  axis.text.x = element_text(size = 15),
  axis.text.y = element_text(size = 15),
  #legend
  legend.title = element_text(size = 20, face = "bold"),
  legend.text = element_text(size = 15)
  )

dev.off()
getwd()
