#randomly assign library-tag combos to PCR products (for 2nd step of 2-step PCR protocol)
#RPK February 2016

#NOTE -- calculate the total number of unique library-tag combinations, and make sure we use this as a maximum for the number of samples we put on the same run.  It's just cleaner than having a particular lib-tag combo used multiple times (e.g., for different loci)
library(googlesheets)
#read in sample names from gsheet
sampleNames=gs_read_csv(gs_title("Samples_OA_good"), ws = "Sequencing_run_Apr_18b")

Ntags=16 #how many tagged primers we have, for the locus of interest
NtotalSamples=nrow(sampleNames)
NusesTag=ceiling(NtotalSamples/Ntags)  #how many times does a given tag get used?
Nloci=1  
Nlibraries=ceiling(NusesTag)*Nloci  #min number of libraries to use

tagNames=paste0("Tag_", seq(1:Ntags))
libNames=paste0("Lib_", LETTERS[1:Nlibraries])

all_combos=expand.grid(tagNames,libNames)  #create all combinations of tags and libraries
all_combos[,3]=paste(all_combos[,1], all_combos[,2], sep="_")  #paste these together into a vector

b=sample(nrow(all_combos))#randomize that vector

c=all_combos[b,]
dim(all_combos)<-c(ceiling(length(all_combos)/Nloci),Nloci)  #reshape the vector into a data frame with 1 column per locus we need for this project



output=data.frame(sampleNames, c[1:NtotalSamples,]) 
	names(output)=c(names(sampleNames), c("Tag", "Lib", "Combo"))

with (output, table(Tag,Lib))	
gs_ws_new(gs_title("Samples_OA_good"), ws_title = "Sequencing_run_Apr_18b", input = output, trim = TRUE)	
write.csv(output, "/Users/rgallego/Google_Drive/Kelly_Lab/Projects/OA_eDNA/Data/Seq_run_Apr_18_with_tags.csv", row.names = F)
output<-read.csv("/Users/rgallego/Google_Drive/Kelly_Lab/Projects/OA_eDNA/Data/Seq_run_Apr_18_with_tags.csv")

output[,"Column"]= gsub("[^[:digit:]]","",output[,2]) #Extract the Clumn from the original plate
output[,"Row"]=as.factor(gsub("[[:digit:]]","",output[,2])) #Extractt the Row from the original plate
table_PCR2=data.frame(matrix(ncol = 12, nrow = 8))          #Empty data.frame with a 96 well plate format
rownames(table_PCR2)=levels(output[,"Row"])
colnames(table_PCR2)=1:12
table_ligation<-table_PCR2
output2<-output
output2$Tag<-as.character(output2$Tag)
output2$Lib<-as.character(output2$Lib)
for (i in LETTERS[1:8]){
  for (j in 1:12){
    well<-paste0(i,j)
    ifelse(well %in% output$Well_in_AMPureXP,
           table_PCR2[i,j]<-output2[output2$Row==i&output2$Column==j,3],table_PCR2[i,j]<-"NA")
    ifelse(well %in% output$Well_in_AMPureXP,
     table_ligation[i,j]<-output2[output2$Row==i&output2$Column==j,4],table_ligation[i,j]<-"NA" )#Find the Lib Column)
  }
}
write.csv(table_PCR2,"~/Google_Drive/Kelly_Lab/Projects/HALO/Data/PCR2.csv", row.names = F)
write.csv(table_ligation,"~/Google_Drive/Kelly_Lab/Projects/HALO/Data/ligation.csv", row.names = F)
##NOw we hac=ve used the updated spreadsheet to add the Qubit concentration of each sample. SO we'll have to load it again
sampleNames=gs_read_csv(gs_title("Samples_OA_good"), ws = "Sequencing_run_Apr_18b")
samplesperlib<- with(sampleNames, rowSums(table(Lib,Tag))) #That tells us how many samples per lib
avgvolpersample<-50/samplesperlib #A way of flagging which samples can give me trouble
ngpersample<-150/samplesperlib#how many ng of DNA per sample into each library

sampleNames<-sampleNames %>% group_by(Lib) %>% mutate(vol_to_pool=round(ngperlib[Lib]/QUBIT_GOOD,1))
#table_for_pooling<-data.frame(matrix(ncol = 12, nrow = 8))   
sampleNames<-sampleNames %>% group_by(Lib) %>% mutate(danger=ifelse(vol_to_pool>avgvolpersample[Lib],"YES","NO"))
sampleNames %>% group_by(Lib) %>% summarise(sum(vol_to_pool))
summary(as.factor(sampleNames$danger))

vol_to_ligation<-table_PCR2
for (i in LETTERS[1:8]){  #This is awful Re think and do properly - maybe a function and then lapply
  for (j in 1:12){
    well<-paste0(i,j)
    vol_to_ligation[i,j]<-""
    ifelse(well %in% sampleNames$Well_in_AMPureXP,
           vol_to_ligation[i,j]<-sampleNames[sampleNames$Row==i & sampleNames$Column==j,"vol_to_pool"],
           vol_to_ligation[i,j]<-"")
    
    ifelse(table_ligation[i,j]=="Lib_G",  #Change across libraries
           print("hi"),vol_to_ligation[i,j]<-"") #Find the Vol to pool column
      

    
  }
}

sampleNames %>% group_by(Lib) %>% summarise(50-sum(vol_to_pool))


rownames(table_PCR2)=levels(output[,"Row"])
colnames(table_PCR2)=1:12
for (i in levels(output[,"Row"])){
  for (j in 1:12){
    table_for_pooling[i,j]<-10*output_df[output$Row==i&output$Column==j,"vol_per_sample_uL"]
    
  }
}