library("dada2")
packageVersion("dada2")


# prepare data information
path <- "/mnt/data/chomic/seqme2020/pl_regenerace/its_znova/dada/"
raw_reads <- sort(list.files("/mnt/data/chomic/seqme2020/pl_regenerace/its_znova/raw", 
                             pattern="*_R1.fq", full.names = TRUE))
sample_names <- sapply(strsplit(basename(raw_reads), "[_]"), `[`, 1)

# inspect quality profiles
plotQualityProfile(raw_reads[1:25])
plotQualityProfile(raw_reads[26:50])
plotQualityProfile(raw_reads[51:80])

plotQualityProfile(filt[1])



# discard reads with unambiguous bases and >2 expected errors and truncate at lenght 245 bp
filt <- file.path(path, "filt", paste0(sample_names, "_filt.fq.gz"))
filt_out <- filterAndTrim(raw_reads, filt, maxN = 0, maxEE = 2, truncLen = 245, multithread = T)

# construct ASV table
err <- learnErrors(filt, multithread = T)
dada <- dada(filt, err=err, multithread=TRUE)
seqtab <- makeSequenceTable(dada, orderBy = "abundance")

# remove chimeras 
seqtab_nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, 
                                        verbose=TRUE)
# write ASVs FASTA file for further taxonomic assignment
write.fasta.dada<-function(dada2, file){
  seqs<-dada2::getSequences(dada2)
  hash<-paste0(">",sapply(seqs, openssl::sha1, USE.NAMES = F))
  write(c(rbind(hash, seqs)),file)
}

write.fasta.dada(seqtab_nochim, paste0(path,"asvs.fa"))

# export ASV table
seqtab_nochim_exp <- seqtab_nochim
colnames(seqtab_nochim_exp) <- sapply(colnames(seqtab_nochim_exp), openssl::sha1, USE.NAMES = F)
rownames(seqtab_nochim_exp) <- paste0("s", gsub("_filt.fq.gz", "", rownames(seqtab_nochim_exp)))
head(seqtab_nochim_exp)[1:6,1:6]

write.table(t(seqtab_nochim), paste0(path,"asv_table.txt"), quote = F)

