library(edgeR)
library(foreach)
library(doParallel)
library(data.table)

variable = commandArgs(TRUE)
counts =  data.frame(fread(paste0(variable[[1]],"/counts.csv"),data.table = FALSE),row.names=1,check.names = FALSE)  
rname_counts = unique(substr(rownames(counts),1,15))
rsums = rowSums(counts)


col_counts = colnames(counts)  # read in counts
# beads = counts[,grep("Beads|BEADS|beads",col_counts)]
beads = counts[,grep("Bead|BEAD|bead",col_counts)] #isolate beads from sample names

beads = beads[,!grepl("high|HIGH|High|Neuro",colnames(beads))] #this only happens for gabes weird stuff

registerDoParallel(12)

##### this calculates the dispersion of the bead only controls
# the below lines are essentially taken straight from the vignette
group <- c(rep(1,dim(beads)[2]))
y <- DGEList(counts=beads, group=group)
keep <- rowSums(cpm(y)>1) >= 2
y <- y[keep, , keep.lib.sizes=FALSE]
y <- calcNormFactors(y)
y <- estimateCommonDisp(y)
common = y$common.dispersion

#####



final = foreach(R = 1:(dim(counts)[2])) %dopar%{
 # final = foreach(R = 1:4) %dopar%{
  print(c(R,"1"))
    if(sum(counts[,R]) == 0){
      counts[,R] = 1
    }
  enrichment = cbind(beads,counts[,R])   # make df of just one sample and bead counts
  group <- c(rep(1,dim(beads)[2]),2) # split into seperate groups
  y <- DGEList(counts=enrichment, group=group) # edgeR 
  y <- calcNormFactors(y)# edgeR 
  y$common.dispersion = common # edgeR 
  # y = estimateTagwiseDisp(y)
  et <- exactTest(y) # edgeR  - this line is what calculates hte acutaly enrichment statistics - it rpeorts both a p -value and a fold change
  print(c(R,"2"))
  et2 = as.data.frame(round(-log10(et$table$PValue),2))
  et2[ et2>250.1] <- 250.1 # sometimes p values get super wonky with super low p -values to we say -log10(x) cant be bigger than 250
  et2 = sign(et$table$logFC)*et2 # set direction of p-value enrichement
  colnames(et2) = col_counts[R]
  rownames(et2) = row.names(et$table)
  
  et3 = as.data.frame(round(2^(et$table$logFC),2)) # same as for et2 but with just log fold change
  colnames(et3) = col_counts[R]
  rownames(et3) = row.names(et$table)
  et4 = list()
  et4[[1]]= et2
  et4[[2]] = et3
  return(et4) #return df of each peptides p-value and log fc
}
print("passed initial")
  final2 = data.frame(matrix(0,nrow = dim(counts)[1],ncol = dim(counts)[2])) # final 2 is matrix of log p values
  colnames(final2) = col_counts
  rownames(final2) = row.names(beads)
  
  final3 = data.frame(matrix(0,nrow = dim(counts)[1],ncol = dim(counts)[2])) # final 3 is matrix of log fc
  colnames(final3) = col_counts
  rownames(final3) = row.names(beads)
  
  for(R in 1:(dim(counts)[2])){
  # for(R in 1:4){
    final2[,R] = final[[R]][[1]]
    final3[,R] = final[[R]][[2]]
    print(paste(R,final2[1,R]))
  }

 
fwrite(as.data.frame(final2),paste0(variable[[2]],"enrichment.csv"),row.names = TRUE)
fwrite(as.data.frame(final3),paste0(variable[[2]],"fold_change.csv"),row.names = TRUE)

file1 = data.frame(fread(paste0(variable[[1]],"/counts.csv")),row.names=1,check.names = FALSE)
file2 = data.frame(fread(paste0(variable[[1]],"/enrichment.csv")),row.names=1,check.names = FALSE)
file3 = data.frame(fread(paste0(variable[[1]],"/fold_change.csv")),row.names=1,check.names = FALSE)
file1[ file1<15 ] <- 0
file1[ file1>=15 ] <- 1
file2[ file2<3 ] <- 0
file2[ file2>=3 ] <- 1
file3[ file3<5 ] <- 0
file3[ file3>=5 ] <- 1
file4 = (file1+file2+file3)/3
file4[ file4<1 ] <- 0

# our final definitation of a hit required the above requiremnts, a counts about 15, log pvalue of 3 and log fc of 5

fwrite(file4,paste0(variable[[2]],"/Hits.csv"),row.names = TRUE)


