library(data.table)
library(stringr)
library(readr)
library(dplyr)
##install necessary R files

####PCA and PRS from  1000 Genomes
pca=read.table('../data/1000G_PCA.txt',header=T)


####Population Model
population_model = glm(ldpred ~ PC1 + PC2 + PC3 + PC4  , data = pca, family = "gaussian")

pca$residual_score = resid(population_model)
pca$residual_score2 = resid(population_model)^2
population_resid_mean = mean(pca$residual_score)
population_resid_sd = sd(pca$residual_score)

####Variance Model
population_var_model <- glm(residual_score2 ~ PC1 + PC2 + PC3 + PC4  , data = pca, family = "gaussian")

generate_adjusted_scores = function(new_data) {
  new_data_adjusted <- new_data %>% mutate(adjusted_score = (ldpred - predict(population_model, new_data))/sqrt(predict(population_var_model, new_data)))
  new_data_adjusted %>% mutate(percentile=pnorm(adjusted_score,0))
}


##############
##SNP Lists for PC projections

map1=fread('../data/1000G_map.txt')
vcf= #input the vcf file

#Extarcting SNP from the single sample imputed files
vvcf=c()

for(i in 1:22){
  map=map1[map1$chromosome==i,]
  tmp=fread(paste0(vcf,i,'.vcf.gz')) ##change as per the input vcf file
  idx=match(paste0(map$chromosome,':',map$physical.pos,'_',map$allele2,'/',map$allele1),paste0(tmp$ID,'_',tmp$REF,'/',tmp$ALT))
  vvcf=rbind(vvcf,tmp[idx,])
}

###Ensure Map1 files and extarcted SNPs are in the order
idx=match(paste0(map1$chromosome,':',map1$physical.pos,'_',map1$allele2,'/',map1$allele1),paste0(vvcf$ID,'_',vvcf$REF,'/',vvcf$ALT))
vvcf=vvcf[idx,]

# RS Find dosage field
tp= str_split_fixed(vvcf$FORMAT, ":", n=4)
tp1= str_split_fixed(vvcf$SAMPLE, ":", n=4)

idx1=which(tp[,1]=='DS')
idx2=which(tp[,2]=='DS')
idx3=which(tp[,3]=='DS')
idx4=which(tp[,4]=='DS')


## RS Dosage of vvcf sites in the order of 1kg map file
dosage=rep(0,nrow(tp))

if(length(idx1)!=0){
  dosage[idx1]=as.numeric(tp1[idx1,1])
}

if(length(idx2)!=0){
  dosage[idx2]=as.numeric(tp1[idx2,2])
}
if(length(idx3)!=0){
  dosage[idx3]=as.numeric(tp1[idx3,3])
}

if(length(idx4)!=0){
  dosage[idx4]=as.numeric(tp1[idx4,4])
}


center=fread('../data/1000G_center.txt')
scale=fread('../data/1000G_scale.txt')
v=fread('../data/1000G_PC1.txt')
d=fread('../data/1000G_lambda.txt')

##Transform Dosage
dosage=2-dosage
#Transform dosage
X9 <- (dosage - center$out.center)/as.numeric(scale$out.scale)

loadings <-  (as.numeric(X9) %*% as.matrix(v))

PC=data.frame(loadings)
for(i in 1:10){
  PC[,i]=PC[,i]/d$out.d[i]
}

#This is the 'raw' PRS from the single sample file, output  of Python PRS script
PC$ldpred=prs_chr   ######PRS from the previous R/Python script
colnames(PC)[1:5]=c('PC1','PC2','PC3','PC4','PC5')
colnames(PC)[6:10]=c('PC6','PC7','PC8','PC9','PC10')

#############PCs Adjusted PRS

new_data <- generate_adjusted_scores(data.frame(PC)) #For the sample
