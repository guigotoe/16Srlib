f1 <- read.table("/home/torres/Documents/Projects/Metagenome/docs/16S/first_info_Femke/FoCus_BSPSPC_Age-Gender_femke_info.txt",header=T,sep="\t",blank.lines.skip=TRUE,na.strings=c("","NA"))
f2 <- read.table("/home/torres/ikmb_storage/Metagenome/16s/2017/4_2017/16s.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.subsample.opti_mcc.unique_list.count.summary",header=F,sep="\t",blank.lines.skip=TRUE,na.strings=c("","NA"))
f3 <- read.table("/home/torres/ikmb_storage/Metagenome/rawdata/16S/MasterFile_plus_LLI03_207.txt",header=T,sep="\t",blank.lines.skip=TRUE,na.strings=c("","NA"))
f3.1 <- subset(f3,readLib==1)[,c(1,3,6,8)]


design <- f2
colnames(design) <- c("ID","LibSize")
CC <- c()
gender <- c()
age <- c()

for (i in design$ID){
  #message(i)
  cc <- as.character(f3.1[f3.1$CoreIDs==i,4])
  if(identical(cc,character(0))) message(i)
  line <- f1[f1$new_id==cc,]
  if (nrow(line)==0) {
    line[1,] <- NA
  }
  CC <- c(CC,cc)
  gender <- c(gender,as.character(line$gender))
  age <- c(age,as.numeric(line$age.years))
}
design$CC <- CC
design$gender <- gender
design$age <- age

techrep <- unlist(lapply(design$ID,function (x) paste(strsplit(as.character(x),'_',fixed=T)[[1]][-1],collapse="_") ))
techrep[techrep==''] <- NA
design$techrep <- techrep
missage <- design[is.na(design$age)|is.na(design$gender),]
write.table(missage,'/home/torres/ikmb_storage/Metagenome/rawdata/16S/missage.txt',quote=F,sep="\t",row.names=F) ## for manual correction using xptools
f4 <- read.table("/home/torres/ikmb_storage/Metagenome/rawdata/16S/missage_filled.txt",header=T,sep="\t",blank.lines.skip=TRUE,na.strings=c("","NA"))
f4$gender <- as.character(f4$gender)
f4$age <- as.numeric(f4$age)
f4$techrep <- as.character(f4$techrep)
f4$CC <- as.character(f4$CC)
for (i in f4$ID) {design[design$ID==i,] <- f4[f4$ID==i,]}

## correcting Central Codes
design$CC <- unlist(lapply(design$CC,function(x) gsub('BBSP','BSP',x,fixed=T) )) 
design$CC <- unlist(lapply(design$CC,function(x) gsub('AGE1','AGE',x,fixed=T) ))

# correcting sample info
design$gender[design$CC=='1165027BSP'] <- 'male' 
design$gender[design$CC=='1165027BSP'] <- 32 

write.table(design,'/home/torres/ikmb_storage/Metagenome/rawdata/16S/design1500.txt',quote=F,sep="\t",row.names=F)



