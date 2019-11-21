#randomly assign library-tag combos to PCR products (for 2nd step of 2-step PCR protocol)
#RPK February 2016

#NOTE -- calculate the total number of unique library-tag combinations, and make sure we use this as a maximum for the number of samples we put on the same run.  It's just cleaner than having a particular lib-tag combo used multiple times (e.g., for different loci)

#read in sample names
sampleNames=read.table("/Users/rpk/Desktop/Willapa_PCR2_samplenames.txt")


Ntags=30 #how many tagged primers we have, for the locus of interest
NtotalSamples=nrow(sampleNames)
NusesTag=ceiling(NtotalSamples/Ntags)  #how many times does a given tag get used?
Nloci=1  
Nlibraries=ceiling(NusesTag)*Nloci  #min number of libraries to use

tagNames=paste0("Tag_", seq(1:Ntags))
libNames=paste0("Lib_", LETTERS[1:Nlibraries])

a=expand.grid(tagNames,libNames)  #create all combinations of tags and libraries
a=paste(a[,1], a[,2], sep="_")  #paste these together into a vector
a=sample(a) #randomize that vector
dim(a)<-c(ceiling(length(a)/Nloci),Nloci)  #reshape the vector into a data frame with 1 column per locus we need for this project



output=data.frame(sampleNames, a[1:NtotalSamples,]) 
	names(output)=c("Sample", c(paste0("Locus", Nloci)))
write.csv(output, "/Users/rpk/GoogleDrive/Kelly_Lab/Projects/AquacultureDiversity/Data/Willapa_PCR2_samplenames_wTagsLibraries.txt")






