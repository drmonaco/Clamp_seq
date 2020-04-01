library(data.table)
library(rBLAST)
library(seqinr)
library(Biostrings)
library(doParallel)
library(foreach)
library(tidyverse)

dbp <- blast("/data/hlarman1/PhIPdb/Studies/Daniel/Phageome/Blast/all_phage.fasta",type="tblastn")
seqp <- readAAStringSet("/data/hlarman1/PhIPdb/Studies/Daniel/Phageome/Blast/PhageomeLa_001_rename.fasta")

#cycle through each patient column by column (goal is so be serialized)

  predp1 <- predict(dbp,seqp[1:2000],BLAST_args = "-evalue .1 -seg no -max_hsps 1 -soft_masking false -word_size 7 -max_target_seqs 100000",custom_format = 
                      "qseqid sseqid pident length mismatch gapopen qstart qend sstart send qseq sseq evalue bitscore nident positive gaps ppos")
 
fwrite(predp1 %>% as.data.frame(),"/data/hlarman1/PhIPdb/Studies/Daniel/Phageome/output_2000.csv")


