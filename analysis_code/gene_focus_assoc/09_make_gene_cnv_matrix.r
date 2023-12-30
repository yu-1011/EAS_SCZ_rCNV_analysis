## read in command line arguments
args <- commandArgs(TRUE)
wdir <- args[1]
input <- args[2]
gene_file <- args[3]
phefile <- args[4] 
region_file <- args[5]
output <- args[6]

print(paste0("Read cnv_file ",input))
print(paste0("Read gene_location ",gene_file))
print(paste0("Read phe_location ",phefile))
print(paste0("Read loci_file ",region_file))

## change to working directory
setwd(wdir)

phe <- read.delim(phefile,h=T,sep=" ")
phe$matchID <- paste(phe$FID,':',phe$IID,sep='')

##get intersection of fam file and phe
fam <- read.table(paste0(input, ".fam"))
colnames(fam) <- c("FID", "IID", "Within-family_ID_of_father", "Within-family_ID_of_mother", "sex", "affect")
fam$matchID <- paste(fam$FID, ':',fam$IID,sep='')
combined <- merge(fam, phe, by='matchID')

## PLINK permutation count (currently ignoring PLINK p-values)
perm <- 100 

## dataset - sort by sizes
dat <- unique(sort(combined$CNV_platform))
data_name <- c('combined',as.character(dat))

##set the start of datasets loop
startpoint=1
if(length(data_name)==2){
	startpoint = 2
}
## write out PLINK ID files 
for (i in 1:1) {
    ID <- combined[combined$CNV_platform==data_name[i],2:3]
    write.table(ID,paste(data_name[i],'.ID',sep=''),col=F,row=F,quo=F,sep='\t')
}

loci <- read.table(region_file)

type_name <- c('allCNV','del','dup')
region_name <- c(as.character(loci$V4))
size_name <- c('>20kb','>100kb','>200kb','>300kb','>400kb','>500kb','>600kb') 
freq_name <- c('Allfreq','singleton','2-5','6-10','11-20','21-40','41-80','81+')

# CNV burden filters
data <- c('',paste(' --keep ', wdir,'/',data_name[2:length(data_name)],'.ID',sep=''))

type <- c('',
          '--cnv-del',
          '--cnv-dup')

region <- c(paste0('--cnv-intersect ', region_file))

size <- c('',
          '--cnv-kb 100',
          '--cnv-kb 200',
          '--cnv-kb 300',
          '--cnv-kb 400',
          '--cnv-kb 500',
          '--cnv-kb 600') 

freq <- c('',
          '--cnv-freq-method2 0.5 --cnv-freq-exclude-above 1',
          '--cnv-freq-method2 0.5 --cnv-freq-exclude-above 5 --cnv-freq-exclude-below 2',
          '--cnv-freq-method2 0.5 --cnv-freq-exclude-above 10 --cnv-freq-exclude-below 6',
          '--cnv-freq-method2 0.5 --cnv-freq-exclude-above 20 --cnv-freq-exclude-below 11',
          '--cnv-freq-method2 0.5 --cnv-freq-exclude-above 40 --cnv-freq-exclude-below 21',
          '--cnv-freq-method2 0.5 --cnv-freq-exclude-above 80 --cnv-freq-exclude-below 41',
          '--cnv-freq-method2 0.5 --cnv-freq-exclude-below 81')


## FILTER sets:
for(a in startpoint:1){
for(b in 1:length(type_name)){
for(c in 1:length(region_name)){
for(d in 1:1){
for(e in 1:1){

## CNV burden loop
##  - .cnv with gene annotation
tmp_data <- list()
X <- 0 ## loop COUNTER
data_set <- NA; region_set <- NA; CNV_type <- NA; CNV_freq <- NA; CNV_size <- NA

## FILTER sets:
X <- X+1
tot <- length(data_name)*length(type_name)*length(region_name)*length(freq_name)*length(size_name)

##Setting up burden loop
## assign specific loop parameters
data2 <- data[a]
type2 <- type[b]
region2 <- region[c]
freq2 <- freq[d]
size2 <- size[e]

## combining CNV type and region filtering
data_set[X] <- data_name[a]
CNV_type[X] <- type_name[b]
region_set[X] <- region_name[c]
CNV_freq[X] <- freq_name[d]
CNV_size[X] <- size_name[e]

## create temporary subdirectory within working directory (will not overwrite existing directory)
final_output_wd <- paste0(output,"/", data_set[X],"_",region_set[X],"_",CNV_type[X],"_",CNV_freq[X],"_",gsub(">","",CNV_size[X]))
print(paste0("output directory: ",final_output_wd))

system(paste0("mkdir -p ",final_output_wd,"/burden_loop"))

## create final subdirectory within output directory (will not overwrite existing directory)

## PART 3.2: Applying PLINK filters


## =============== Filter by data, CNV type, and CNV regions first

  system(paste("/stanley/huang_lab/home/ychen/software/plink-1.07-x86_64/plink --noweb --cfile ",input," ",data2," ",type2," ",region2," --cnv-write --out ",final_output_wd,"/burden_loop/type_loop",sep=""))
  
  ## gene cnv-make-map
  system(paste("/stanley/huang_lab/home/ychen/software/plink-1.07-x86_64/plink --noweb --cfile ",final_output_wd,"/burden_loop/type_loop --cnv-make-map --out ",final_output_wd,"/burden_loop/type_loop",sep=""))


  ## --- Frequency pruning    
  system(paste("/stanley/huang_lab/home/ychen/software/plink-1.07-x86_64/plink --noweb --cfile ",final_output_wd,"/burden_loop/type_loop ",freq2," --cnv-write --out ",final_output_wd,"/burden_loop/freq_loop",sep=""))
  ## gene cnv-make-map
  system(paste("/stanley/huang_lab/home/ychen/software/plink-1.07-x86_64/plink --noweb --cfile ",final_output_wd,"/burden_loop/freq_loop --cnv-make-map --out ",final_output_wd,"/burden_loop/freq_loop",sep=""))

  ## --- Size pruning
  system(paste("/stanley/huang_lab/home/ychen/software/plink-1.07-x86_64/plink --noweb --cfile  ",final_output_wd,"/burden_loop/freq_loop ",size2," --cnv-write --out ",final_output_wd,"/burden_loop/size_loop",sep=""))
  system(paste("/stanley/huang_lab/home/ychen/software/plink-1.07-x86_64/plink --noweb --cfile  ",final_output_wd,"/burden_loop/size_loop --cnv-make-map --out  ",final_output_wd,"/burden_loop/size_loop",sep=""))
 
  ## Adding gene/exon count separately
  system(paste("/stanley/huang_lab/home/ychen/software/plink-1.07-x86_64/plink --noweb --cfile  ",final_output_wd,"/burden_loop/size_loop --cnv-indiv-perm --mperm ",perm," --cnv-count ",gene_file," --out  ",final_output_wd,"/burden_loop/burden_loop",sep=""))
  
  ## summary final cnv dataset
  system(paste("/stanley/huang_lab/home/ychen/software/plink-1.07-x86_64/plink --noweb --cfile  ",final_output_wd,"/burden_loop/size_loop --cnv-freq-method2 0.5 --cnv-write-freq --cnv-track "," --cnv-count ",gene_file," --out ",final_output_wd,"/",region_set[X], sep=""))
  system(paste("/stanley/huang_lab/home/ychen/software/plink-1.07-x86_64/plink --noweb --cfile  ",final_output_wd,"/burden_loop/size_loop --cnv-freq-method2 0.5 --cnv-write-freq --cnv-write "," --cnv-count ",gene_file," --out ",final_output_wd,"/",region_set[X], sep=""))


  ## remove the temporary subdirectory 
  system(paste0("rm -r  ",final_output_wd,"/burden_loop"))
  system(paste0("rm  ",final_output_wd,"/*.fam"))
}
}
}
}
}

