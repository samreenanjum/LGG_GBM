## separate the protected files based on SS Code and LOH condition
## Script to call the function at the end

split.SS = function(input.filename, file.num){
  title.line <- readLines(paste0(path, "broad.mit.edu_LGG.IlluminaGA_DNASeq_Cont_curated.Level_2.1.2.0/", input.filename), n=33)  
  sample.name = input.filename
  input.table = read.table(paste0(path, "broad.mit.edu_LGG.IlluminaGA_DNASeq_Cont_curated.Level_2.1.2.0/", input.filename), comment.char = "#", stringsAsFactors = F)
  colnames(input.table) = c("chrom", "pos", "id", "ref", "alt", "qual", "filter", "info", "format", "normal", "tumor")
  
  ##Split into somatic and germline
  somatic.calls = input.table[which(grepl(x = input.table$info, pattern = ";SS=Somatic;")),]
  germline.calls = input.table[which(grepl(x = input.table$info, pattern =  ";SS=Germline;")),]
  
  pos = strsplit(somatic.calls$format[1], split = ":")
  pos.freq = which(pos[[1]]=="FA");   pos.dp = which(pos[[1]]=="DP");  
  
  ## Germline Filter
  Germline.TumorDP = sapply(strsplit(germline.calls$tumor, split = ":"), "[[", pos.dp)
  Germline.NormalDP = as.numeric(sapply(strsplit(germline.calls$normal, split = ":"), "[[", pos.dp))
  germline.reqd.ind = which(Germline.NormalDP>=60)# & TumorVarFreq >=10)
  germline.calls.final = germline.calls[germline.reqd.ind,]
  
  ## Extract details = VAF and DP
  TumorVarFreq = sapply(strsplit(somatic.calls$tumor, split = ":"), "[[", pos.freq)
  NormalVarFreq = sapply(strsplit(somatic.calls$normal, split = ":"), "[[", pos.freq)
  
  TumorDP = sapply(strsplit(somatic.calls$tumor, split = ":"), "[[", pos.dp)
  NormalDP = as.numeric(sapply(strsplit(somatic.calls$normal, split = ":"), "[[", pos.dp))
  
  TumorVarFreq = as.double(TumorVarFreq)*100
  NormalVarFreq = as.double(NormalVarFreq)*100
  TumorVarFreq[which(is.na(TumorVarFreq))] = 0
  
  ## Somatic Filter =  normal DP >=60 
  reqd.ind = which(NormalDP>=60)#
  somatic.calls.filt = somatic.calls[reqd.ind,]
  TumorVarFreq.filt = TumorVarFreq[reqd.ind]
  NormalVarFreq.filt = NormalVarFreq[reqd.ind]
  
  ## Extract LOH and mark remaining as Somatic
  loh.ind = which(abs(TumorVarFreq.filt - 50) > abs(NormalVarFreq.filt-50))
  LOH.calls = somatic.calls.filt[loh.ind,]
  
  somatic.calls.final = somatic.calls.filt[-loh.ind,]

  ## Output files  
  output.g.filename = paste0(path,"Broad.Filtered.Reads/Germline/", sample.name, ".germline.vcf")
  output.s.filename = paste0(path,"Broad.Filtered.Reads/Somatic/", sample.name, ".somatic.vcf")
  output.l.filename = paste0(path,"Broad.Filtered.Reads/LOH/", sample.name, ".LOH.vcf")
  
  write(title.line, file = output.g.filename)
  write.table(germline.calls.final, output.g.filename, row.names=F,  quote = F , col.names = F, append = T, sep = "\t")
  
  write(title.line, file = output.s.filename)
  write.table(somatic.calls.final, output.s.filename, row.names=F,  quote = F , col.names = F, append = T, sep = "\t")
  
  write(title.line, file = output.l.filename)
  write.table(LOH.calls, output.l.filename, row.names=F,  quote = F , col.names = F, append = T, sep = "\t")
   
}

path = "~/LGG-GBM/Protected/Broad/"
files = list.files(path = paste0(path, "broad.mit.edu_LGG.IlluminaGA_DNASeq_Cont_curated.Level_2.1.2.0/"), pattern = ".vcf")

for (file.num in 1:length(files)){
  print(file.num)
  input.filename = files[file.num]#"TCGA-S9-A6TW-01A-12D-A32B-08.somatic.snp.vcf"
  split.SS(input.filename, file.num)
}

