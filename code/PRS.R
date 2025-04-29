library(data.table)
library(readr)
library(stringr)

prs=c()
prs_chr=0

vcf='genome.imputed.chr' ###Change the vcf file as per the input
file_path=fread('../data/trans_prs_Nov_19.txt')

for(i in 1:22){
  file=file_path[file_path$chr==i,]
  tmp=fread(paste0(vcf,i,'.vcf.gz'))  ###change the vcf file
  idx=match(file$newid,paste0(tmp$ID,':',tmp$REF,':',tmp$ALT))
  #prs=rbind(prs,tmp[idx,])

  tpp=tmp[idx,]
  tp= str_split_fixed(tpp$FORMAT, ":", n=4)
  tp1= str_split_fixed(tpp$SAMPLE, ":", n=4)

  idx1=which(tp[,1]=='DS')
  idx2=which(tp[,2]=='DS')
  idx3=which(tp[,3]=='DS')
  idx4=which(tp[,4]=='DS')


  fprs=rep(0,nrow(tp))

  if(length(idx1)!=0){
    fprs[idx1]=as.numeric(tp1[idx1,1])
  }

  if(length(idx2)!=0){
    fprs[idx2]=as.numeric(tp1[idx2,2])
  }
  if(length(idx3)!=0){
    fprs[idx3]=as.numeric(tp1[idx3,3])
  }

  if(length(idx4)!=0){
    fprs[idx4]=as.numeric(tp1[idx4,4])
  }

  #idx=match(paste0(tprs$ID,':',tprs$REF,':',tprs$ALT),file_path$newid)
  final_prs1 <-  (as.numeric(fprs) %*% as.matrix(file$beta_grid4))


  prs_chr=prs_chr+final_prs1

}
