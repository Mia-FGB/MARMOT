cat("Script execution started.\n")

#Script mainly from https://github.com/mgiolai/crop_airseq/tree/main/prepare-PHIbase-genomes 
library(tidyverse)
library(taxizedb)
library(RCurl)
library(taxonomizr)


#prepareDatabase('accessionTaxa.sql') - #only need to do this once
taxDB <- "accessionTaxa.sql"

# DEFRA database https://planthealthportal.defra.gov.uk/pests-and-diseases/uk-plant-health-risk-register/downloadEntireRiskRegister.cfm----------------------------------------------------------
#Download the most up to date version each time

#-----Importing the data, this path needs to change each time
cat("Loading data from Risk_Register and PHIbase...\n")
defra <- read.csv("~/Library/CloudStorage/OneDrive-NorwichBioScienceInstitutes/Pathogen_Database/Pathogen_Database_080524/Risk_Register_05_24.csv")

#------Don't want insect/nematode/plant/mite pests in this database
remove <-  c("Insect", "Mite", "Nematode", "Plant")
defra <- defra %>% 
  filter(!(Type.of.pest %in% remove))

#----Printing a list of the species names 
defra_pest_name <- defra %>% 
  select(Pest.Name, Type.of.pest)

#----Option to add in a row for additional pathogens of interest: 
defra_pest_name[nrow(defra_pest_name) + 1,] = c('Erysiphe necator', 'Fungus')

#----Get TaxaID from species name
defra_name_ID <- name2taxid(defra_pest_name$Pest.Name, db = "ncbi", verbose = TRUE, out_type ="summary")
defraID <- defra_name_ID %>% 
  select(id)

# PHIbase https://github.com/PHI-base/data/tree/master/releases -----------------------------------------------------------------
#Download the most up to date versio each time 

#-----Importing the data, this path needs to change each time
dfPhibase <- read.delim("~/Library/CloudStorage/OneDrive-NorwichBioScienceInstitutes/Pathogen_Database/Pathogen_Database_080524/phi-base_0524.csv",sep=",", fill = T)

#---------Get PHIbase species ids 
dfPhibase <- dfPhibase %>% 
  select(Pathogen_NCBI_species_Taxonomy.ID) %>% 
  rename(id = Pathogen_NCBI_species_Taxonomy.ID )

dfPhibase$id <- as.numeric(dfPhibase$id)


#Combine the two list of species & taxaIDs together----------------
cat("Combining data from PHIbase and DEFRA databases...\n")
combdf <- rbind(dfPhibase, defraID)

#Getting unique taxaIDs, and remove na values
un_combdf <- unique(combdf) %>% 
  na.omit(combdf$id)

#Taxonomise with taxonomizr phylogeny
cat("Taxonomising analysis...\n")
dfTax_db <- data.frame(getTaxonomy(unique(un_combdf$id),taxDB)) #get lineage from taxaID
dfTax_db$taxid <- as.numeric(gsub(" ","",row.names(dfTax_db))) #taking taxa id to own column from row name
df_name_ID <- merge(un_combdf,dfTax_db,by.x="id",by.y="taxid",all=T) #add taxonomy to db
write.csv(df_name_ID,"~/Library/CloudStorage/OneDrive-NorwichBioScienceInstitutes/Pathogen_Database/Pathogen_Database_080524/Pathogen_Database_080524.taxonomised.csv")

# Download Refseq & Genbank tables ----------------------------------------

#Download refseq and genbank tables - Need to do this each time so they are up to date
download.file("https://ftp.ncbi.nlm.nih.gov/genomes/ASSEMBLY_REPORTS/assembly_summary_refseq.txt","assembly_summary_refseq.txt")
download.file("https://ftp.ncbi.nlm.nih.gov/genomes/ASSEMBLY_REPORTS/assembly_summary_genbank.txt","assembly_summary_genbank.txt")

#Read in refseq dataframe
refseq <- read.delim("assembly_summary_refseq.txt",header=F,sep="\t",skip=2,quote="")
colnames(refseq) <- c("assembly_accession","bioproject","biosample","wgs_master","refseq_category","taxid","species_taxid","organism_name","infraspecific_name","isolate","version_status","assembly_level","release_type","genome_rep","seq_rel_date","asm_name","submitter","gbrs_paired_asm","paired_asm_comp","ftp_path","excluded_from_refseq","relation_to_type_material","asm_not_live_date")
# Remove columns with NA names
refseq <- refseq[!is.na(names(refseq))]

#Read in genbank dataframe
genbank <- read.delim("assembly_summary_genbank.txt",header=F,sep="\t",skip=2,quote="")
colnames(genbank) <- c("assembly_accession","bioproject","biosample","wgs_master","refseq_category","taxid","species_taxid","organism_name","infraspecific_name","isolate","version_status","assembly_level","release_type","genome_rep","seq_rel_date","asm_name","submitter","gbrs_paired_asm","paired_asm_comp","ftp_path","excluded_from_refseq","relation_to_type_material","asm_not_live_date")
# Remove columns with NA names
genbank <- genbank[!is.na(names(genbank))]


#Remove entries from genbank that are present in refseq to avoid redundancies
genbank <- genbank %>% dplyr::filter(!assembly_accession %in% refseq$gbrs_paired_asm)

#Merge the assembly databases
ref_gen <- rbind(refseq,genbank)

#Taxonomise with taxonomizr phylogeny
dfTax_ref_gen <- data.frame(getTaxonomy(unique(ref_gen$species_taxid),taxDB))
dfTax_ref_gen$taxid <- as.numeric(gsub(" ","",row.names(dfTax_ref_gen)))
ref_gen <- merge(ref_gen,dfTax_ref_gen,by.x="species_taxid",by.y="taxid",all=T)


##FILTER THE TABLES----------------

#Filter for species present among the PHIbase and DEFRA pathogens
ref_gen_match <- ref_gen %>%
  dplyr::filter(species %in% c(df_name_ID$species))

#Only take rows with an https entry
ref_gen_match <- ref_gen_match[grep("https://",ref_gen_match$ftp_path),]

#Export the filtered and taxonomised refseq genbank dataframe
write.table(ref_gen_match,"~/Library/CloudStorage/OneDrive-NorwichBioScienceInstitutes/Pathogen_Database/Pathogen_Database_080524/refseq.genbank.taxonomised.080524.txt",sep="\t",row.names=F,quote=F)


#Loop 1 - Reference genome
#Query if there is exactly one reference genome, if not, determine the longest
ref_gen_matchAccessions <- c()

cat("Query if there is exactly one reference genome, if not, determine the longest \n")
for(sp in setdiff(sort(unique(ref_gen_match$species)),ref_gen_matchAccessions$species)){
  print(sp)
  #Subset for species and reference genomes
  x <- ref_gen_match %>% dplyr::filter(species==sp) %>% dplyr::filter(refseq_category=="reference genome")
  
  #Select the newest assembly version of each accession
  x$v <- as.numeric(gsub(".","",gsub("GCF_","",gsub("GCA_","",x$assembly_accession)),fixed=T))
  x$acc <- sapply(strsplit(as.character(x$assembly_accession),".",fixed=T),'[[',1)
  x$acc <- gsub("GCF","",gsub("GCA","",x$acc))
  x <- x %>%
    dplyr::group_by(acc) %>%
    dplyr::filter(v == max(v)) %>%
    dplyr::ungroup()
  
  #If there is only one accession for this species, add the accession to the list
  if(length(unique(x$assembly_accession))==1){
    ref_gen_matchAccessions <- rbind(ref_gen_matchAccessions,
                          data.frame("accession"=x$assembly_accession,"species"=sp,"type"="Reference Genome"))
    #If there is more than one accession for this species, determine the largest genome
  }else if(length(unique(x$assembly_accession))>1){
    print("Querying genome sizes")
    tmp <- data.frame()
    for(url in x$ftp_path){
      acc <- x %>% dplyr::filter(ftp_path==url) %>% dplyr::pull(assembly_accession)
      url <- paste0(url,"/",tail(strsplit(url,"/",fixed=T)[[1]],n=1),"_assembly_stats.txt")
      download.file(url,"tmp.txt",quiet=T)
      size <- as.numeric(tail(strsplit(grep("all\tall\tall\tall\ttotal-length",readLines("tmp.txt"),value=T),"\t",fixed=T)[[1]],n=1))
      tmp <- rbind(tmp,
                   data.frame("accession"=acc,"size"=size))
      unlink("tmp.txt")
    }
    #Determine the largest accession
    ref_gen_matchAccessions <- rbind(ref_gen_matchAccessions,
                          data.frame("accession"=tmp %>% dplyr::filter(size==max(size)) %>% dplyr::pull(accession),"species"=sp,"type"="Reference Genome"))
  }
}


#Loop 2 - Representative Genome
#Query if there is exactly one representative genome, if not, determine the longest
cat("Query if there is exactly one representative genome, if not, determine the longest \n")

for(sp in setdiff(sort(unique(ref_gen_match$species)),ref_gen_matchAccessions$species)){
  print(sp)
  #Subset for species and representative genomes
  x <- ref_gen_match %>% dplyr::filter(species==sp) %>% dplyr::filter(refseq_category=="representative genome")
  
  #Select the newest assembly version of each accession
  x$v <- as.numeric(gsub(".","",gsub("GCF_","",gsub("GCA_","",x$assembly_accession)),fixed=T))
  x$acc <- sapply(strsplit(as.character(x$assembly_accession),".",fixed=T),'[[',1)
  x$acc <- gsub("GCF","",gsub("GCA","",x$acc))
  x <- x %>%
    dplyr::group_by(acc) %>%
    dplyr::filter(v == max(v)) %>%
    dplyr::ungroup()
  
  #If there is only one accession for this species, add the accession to the list
  if(length(unique(x$assembly_accession))==1){
    ref_gen_matchAccessions <- rbind(ref_gen_matchAccessions,
                          data.frame("accession"=x$assembly_accession,"species"=sp,"type"="Representative Genome"))
    #If there is more than one accession for this species, determine the largest genome
  }else if(length(unique(x$assembly_accession))>1){
    print("Querying genome sizes")
    tmp <- data.frame()
    for(url in x$ftp_path){
      acc <- x %>% dplyr::filter(ftp_path==url) %>% dplyr::pull(assembly_accession)
      url <- paste0(url,"/",tail(strsplit(url,"/",fixed=T)[[1]],n=1),"_assembly_stats.txt")
      download.file(url,"tmp.txt",quiet=T)
      size <- as.numeric(tail(strsplit(grep("all\tall\tall\tall\ttotal-length",readLines("tmp.txt"),value=T),"\t",fixed=T)[[1]],n=1))
      tmp <- rbind(tmp,
                   data.frame("accession"=acc,"size"=size))
      unlink("tmp.txt")
    }
    #Determine the largest accession
    ref_gen_matchAccessions <- rbind(ref_gen_matchAccessions,
                          data.frame("accession"=tmp %>% dplyr::filter(size==max(size)) %>% dplyr::pull(accession),"species"=sp,"type"="Representative Genome"))
  }
}

#Loop 3 - Complete Genome
#Query if there is exactly one complete genome, if not, determine the longest
cat("Query if there is exactly one complete genome, if not, determine the longest \n")

for(sp in setdiff(sort(unique(ref_gen_match$species)),ref_gen_matchAccessions$species)){
  print(sp)
  #Subset for species and complete genomes
  x <- ref_gen_match %>% dplyr::filter(species==sp) %>% dplyr::filter(assembly_level=="Complete Genome")
  
  #Select the newest assembly version of each accession
  x$v <- as.numeric(gsub(".","",gsub("GCF_","",gsub("GCA_","",x$assembly_accession)),fixed=T))
  x$acc <- sapply(strsplit(as.character(x$assembly_accession),".",fixed=T),'[[',1)
  x$acc <- gsub("GCF","",gsub("GCA","",x$acc))
  x <- x %>%
    dplyr::group_by(acc) %>%
    dplyr::filter(v == max(v)) %>%
    dplyr::ungroup()
  
  #If there is only one accession for this species, add the accession to the list
  if(length(unique(x$assembly_accession))==1){
    ref_gen_matchAccessions <- rbind(ref_gen_matchAccessions,
                          data.frame("accession"=x$assembly_accession,"species"=sp,"type"="Complete Genome"))
    #If there is more than one accession for this species, determine the largest genome
  }else if(length(unique(x$assembly_accession))>1){
    print("Querying genome sizes")
    tmp <- data.frame()
    for(url in x$ftp_path){
      acc <- x %>% dplyr::filter(ftp_path==url) %>% dplyr::pull(assembly_accession)
      url <- paste0(url,"/",tail(strsplit(url,"/",fixed=T)[[1]],n=1),"_assembly_stats.txt")
      download.file(url,"tmp.txt",quiet=T)
      size <- as.numeric(tail(strsplit(grep("all\tall\tall\tall\ttotal-length",readLines("tmp.txt"),value=T),"\t",fixed=T)[[1]],n=1))
      tmp <- rbind(tmp,
                   data.frame("accession"=acc,"size"=size))
      unlink("tmp.txt")
    }
    #Determine the largest accession
    ref_gen_matchAccessions <- rbind(ref_gen_matchAccessions,
                          data.frame("accession"=tmp %>% dplyr::filter(size==max(size)) %>% dplyr::pull(accession),"species"=sp,"type"="Complete Genome"))
  }
}

#Loop 4 - Longest genome
#Query if there is exactly one genome, if not, determine the longest
cat("Query if there is exactly one genome, if not, determine the longest \n")

for(sp in setdiff(sort(unique(ref_gen_match$species)),ref_gen_matchAccessions$species)){
  print(sp)
  #Subset for species
  x <- ref_gen_match %>% dplyr::filter(species==sp)
  
  #Select the newest assembly version of each accession
  x$v <- as.numeric(gsub(".","",gsub("GCF_","",gsub("GCA_","",x$assembly_accession)),fixed=T))
  x$acc <- sapply(strsplit(as.character(x$assembly_accession),".",fixed=T),'[[',1)
  x$acc <- gsub("GCF","",gsub("GCA","",x$acc))
  x <- x %>%
    dplyr::group_by(acc) %>%
    dplyr::filter(v == max(v)) %>%
    dplyr::ungroup()
  
  #If there is only one accession for this species, add the accession to the list
  if(length(unique(x$assembly_accession))==1){
    ref_gen_matchAccessions <- rbind(ref_gen_matchAccessions,
                          data.frame("accession"=x$assembly_accession,"species"=sp,"type"="Other"))
    #If there is more than one accession for this species, determine the largest genome
  }else if(length(unique(x$assembly_accession))>1){
    print("Querying genome sizes")
    tmp <- data.frame()
    for(url in x$ftp_path){
      acc <- x %>% dplyr::filter(ftp_path==url) %>% dplyr::pull(assembly_accession)
      url <- paste0(url,"/",tail(strsplit(url,"/",fixed=T)[[1]],n=1),"_assembly_stats.txt")
      download.file(url,"tmp.txt",quiet=T)
      size <- as.numeric(tail(strsplit(grep("all\tall\tall\tall\ttotal-length",readLines("tmp.txt"),value=T),"\t",fixed=T)[[1]],n=1))
      tmp <- rbind(tmp,
                   data.frame("accession"=acc,"size"=size))
      unlink("tmp.txt")
    }
    #Determine the largest accession
    ref_gen_matchAccessions <- rbind(ref_gen_matchAccessions,
                          data.frame("accession"=tmp %>% dplyr::filter(size==max(size)) %>% dplyr::pull(accession),"species"=sp,"type"="Other"))
  }
}


write.csv(ref_gen_matchAccessions,"~/Library/CloudStorage/OneDrive-NorwichBioScienceInstitutes/Pathogen_Database/Pathogen_Database_080524/downloaded_accessions_08054.csv")

# Download file -----------------------------------------------------------

ref_gen_match <- ref_gen_match %>% dplyr::filter(assembly_accession %in% ref_gen_matchAccessions$accession)
write.csv(ref_gen_match,"~/Library/CloudStorage/OneDrive-NorwichBioScienceInstitutes/Pathogen_Database/Pathogen_Database_080524/downloaded_refseq_genbank_080524.csv")

#Add download links
ref_gen_match$dlLink <- paste0(ref_gen_match$ftp_path,"/",paste0(sapply(strsplit(ref_gen_match$ftp_path,"/",fixed=T),tail,1),"_genomic.fna.gz"))
ref_gen_match$dlLinkMD5 <- gsub("ftp:","https:",paste0(ref_gen_match$ftp_path,"/md5checksums.txt"))

#Export download file
write.table(ref_gen_match %>% dplyr::select(species,species_taxid,assembly_accession,dlLinkMD5,dlLink),file="~/Library/CloudStorage/OneDrive-NorwichBioScienceInstitutes/Pathogen_Database/download.txt",sep="\t",row.names=F,col.names=F,quote=F)
cat("Script execution completed.\n")