library(stringr)
library(tidyverse)



######### lipid peroxidation #########
### Non-alcoholic fatty liver disease ###
#read file as data frame
NAF_liver_disease <- read.delim2(file = "Non-alcoholic fatty liver disease.txt", header = FALSE, sep = "\n",dec = ",")
# read as character vector
NAF_liver_disease <- sapply(NAF_liver_disease, as.character) 

#regex for gene names
gene_names <- str_match(NAF_liver_disease, "\\(RefSeq\\)\\s*(.*?)\\s*,")
#remove NA's
gene_names <- na.omit(gene_names) 

# select column with gene names
gene_names <- as.data.frame(gene_names) #needs to be a df otherwise function select will not work
gene_names<-gene_names %>% select(V2)

#write the gene list to csv
write.csv(x=gene_names, file="NAF_liver_disease_gene_names")

### Ferroptosis ###
#read file as data frame
ferroptosis <- read.delim2(file = "ferroptosis.txt", header = FALSE, sep = "\n",dec = ",")
# read as character vector
ferroptosis <- sapply(ferroptosis, as.character) 

#regex for gene names
gene_names <- str_match(ferroptosis, "\\(RefSeq\\)\\s*(.*?)\\s*,")
#remove NA's
gene_names <- na.omit(gene_names) 

# select column with gene names
gene_names <- as.data.frame(gene_names) #needs to be a df otherwise function select will not work
gene_names<-gene_names %>% select(V2)

#write the gene list to csv
write.csv(x=gene_names, file="ferroptosis_gene_names")

### Alcoholic liver disease ###
#read file as data frame
alc_liver_disease <- read.delim2(file = "Alcoholic liver disease.txt", header = FALSE, sep = "\n",dec = ",")
# read as character vector
alc_liver_disease <- sapply(alc_liver_disease, as.character) 

#regex for gene names
gene_names <- str_match(alc_liver_disease, "\\(RefSeq\\)\\s*(.*?)\\s*,")
#remove NA's
gene_names <- na.omit(gene_names) 

# select column with gene names
gene_names <- as.data.frame(gene_names) #needs to be a df otherwise function select will not work
gene_names<-gene_names %>% select(V2)

#write the gene list to csv
write.csv(x=gene_names, file="alc_liver_disease_gene_names")

### Alzheimer disease ###
#read file as data frame
alzheimer_disease <- read.delim2(file = "Alzheimer disease.txt", header = FALSE, sep = "\n",dec = ",")
# read as character vector
alzheimer_disease <- sapply(alzheimer_disease, as.character) 

#regex for gene names
gene_names <- str_match(alzheimer_disease, "\\(RefSeq\\)\\s*(.*?)\\s*,")
#remove NA's
gene_names <- na.omit(gene_names) 

# select column with gene names
gene_names <- as.data.frame(gene_names) #needs to be a df otherwise function select will not work
gene_names<-gene_names %>% select(V2)

#write the gene list to csv
write.csv(x=gene_names, file="alzheimer_disease_gene_names")

### Chemical carcinogenesis - DNA adducts ###
#read file as data frame
chem_carcinogenesis <- read.delim2(file = "Chemical carcinogenesis - DNA adducts.txt", header = FALSE, sep = "\n",dec = ",")
# read as character vector
chem_carcinogenesis <- sapply(chem_carcinogenesis, as.character) 

#regex for gene names
gene_names <- str_match(chem_carcinogenesis, "\\(RefSeq\\)\\s*(.*?)\\s*,")
#remove NA's
gene_names <- na.omit(gene_names) 

# select column with gene names
gene_names <- as.data.frame(gene_names) #needs to be a df otherwise function select will not work
gene_names<-gene_names %>% select(V2)

#write the gene list to csv
write.csv(x=gene_names, file="chem_carcinogenesis_gene_names")



######### redox metabolism #########
### glutathion metabolism ###
#read file as data frame
glutathion_metabolism <- read.delim2(file = "glutathion_metabolism.txt", header = FALSE, sep = "\n",dec = ",")
# read as character vector
glutathion_metabolism <- sapply(glutathion_metabolism, as.character) 

#regex for gene names
gene_names <- str_match(glutathion_metabolism, "\\(RefSeq\\)\\s*(.*?)\\s*,")
#remove NA's
gene_names <- na.omit(gene_names) 

# select column with gene names
gene_names <- as.data.frame(gene_names) #needs to be a df otherwise function select will not work
gene_names<-gene_names %>% select(V2)

#write the gene list to csv
write.csv(x=gene_names, file="glutathion_metabolism_gene_names")

### peroxisome ###
#read file as data frame
peroxisome <- read.delim2(file = "peroxisome.txt", header = FALSE, sep = "\n",dec = ",")
# read as character vector
peroxisome <- sapply(peroxisome, as.character) 

#regex for gene names
gene_names <- str_match(peroxisome, "\\(RefSeq\\)\\s*(.*?)\\s*,")
#remove NA's
gene_names <- na.omit(gene_names) 

# select column with gene names
gene_names <- as.data.frame(gene_names) #needs to be a df otherwise function select will not work
gene_names<-gene_names %>% select(V2)

#write the gene list to csv
write.csv(x=gene_names, file="peroxisome_gene_names")

### chemical carcinogenesis ROS ###
#read file as data frame
chemical_carcinogenesis_ROS <- read.delim2(file = "chemical_carcinogenesis_ROS.txt", header = FALSE, sep = "\n",dec = ",")
# read as character vector
chemical_carcinogenesis_ROS <- sapply(chemical_carcinogenesis_ROS, as.character) 

#regex for gene names
gene_names <- str_match(chemical_carcinogenesis_ROS, "\\(RefSeq\\)\\s*(.*?)\\s*,")
#remove NA's
gene_names <- na.omit(gene_names) 

# select column with gene names
gene_names <- as.data.frame(gene_names) #needs to be a df otherwise function select will not work
gene_names<-gene_names %>% select(V2)

#write the gene list to csv
write.csv(x=gene_names, file="chemical_carcinogenesis_ROS_gene_names")

### mineral absorption ###
#read file as data frame
mineral_absorption <- read.delim2(file = "mineral_absorption.txt", header = FALSE, sep = "\n",dec = ",")
# read as character vector
mineral_absorption <- sapply(mineral_absorption, as.character) 

#regex for gene names
gene_names <- str_match(mineral_absorption, "\\(RefSeq\\)\\s*(.*?)\\s*,")
#remove NA's
gene_names <- na.omit(gene_names) 

# select column with gene names
gene_names <- as.data.frame(gene_names) #needs to be a df otherwise function select will not work
gene_names<-gene_names %>% select(V2)

#write the gene list to csv
write.csv(x=gene_names, file="mineral_absorption_gene_names")

### fatty acid degradation ###
#read file as data frame
fatty_acid_degradation <- read.delim2(file = "fatty_acid_degradation.txt", header = FALSE, sep = "\n",dec = ",")
# read as character vector
fatty_acid_degradation <- sapply(fatty_acid_degradation, as.character) 

#regex for gene names
gene_names <- str_match(fatty_acid_degradation, "\\(RefSeq\\)\\s*(.*?)\\s*,")
#remove NA's
gene_names <- na.omit(gene_names) 

# select column with gene names
gene_names <- as.data.frame(gene_names) #needs to be a df otherwise function select will not work
gene_names<-gene_names %>% select(V2)

#write the gene list to csv
write.csv(x=gene_names, file="fatty_acid_degradation_gene_names")

### protein processing in ER ###
#read file as data frame
protein_processing_in_ER <- read.delim2(file = "protein_processing_in_ER.txt", header = FALSE, sep = "\n",dec = ",")
# read as character vector
protein_processing_in_ER <- sapply(protein_processing_in_ER, as.character) 

#regex for gene names
gene_names <- str_match(protein_processing_in_ER, "\\(RefSeq\\)\\s*(.*?)\\s*,")
#remove NA's
gene_names <- na.omit(gene_names) 

# select column with gene names
gene_names <- as.data.frame(gene_names) #needs to be a df otherwise function select will not work
gene_names<-gene_names %>% select(V2)

#write the gene list to csv
write.csv(x=gene_names, file="protein_processing_in_ER_gene_names")



######### iron metabolism #########
### porphyrin metabolism ###
#read file as data frame
porphyrin_metabolism <- read.delim2(file = "porphyrin_metabolism.txt", header = FALSE, sep = "\n",dec = ",")
# read as character vector
porphyrin_metabolism <- sapply(porphyrin_metabolism, as.character) 

#regex for gene names
gene_names <- str_match(porphyrin_metabolism, "\\(RefSeq\\)\\s*(.*?)\\s*,")
#remove NA's
gene_names <- na.omit(gene_names) 

# select column with gene names
gene_names <- as.data.frame(gene_names) #needs to be a df otherwise function select will not work
gene_names<-gene_names %>% select(V2)

#write the gene list to csv
write.csv(x=gene_names, file="porphyrin_metabolism_gene_names")

### necroptosis ###
necroptosis <- read.delim2(file = "necroptosis.txt", header = FALSE, sep = "\n",dec = ",")
# read as character vector
necroptosis <- sapply(necroptosis, as.character) 

#regex for gene names
gene_names <- str_match(necroptosis, "\\(RefSeq\\)\\s*(.*?)\\s*,")
#remove NA's
gene_names <- na.omit(gene_names) 

# select column with gene names
gene_names <- as.data.frame(gene_names) #needs to be a df otherwise function select will not work
gene_names<-gene_names %>% select(V2)

#write the gene list to csv
write.csv(x=gene_names, file="necroptosis_gene_names")

### TGF_beta_signaling_pathway ###
TGF_beta_signaling_pathway <- read.delim2(file = "TGF-beta_signaling_pathway.txt", header = FALSE, sep = "\n",dec = ",")
# read as character vector
TGF_beta_signaling_pathway <- sapply(TGF_beta_signaling_pathway, as.character) 

#regex for gene names
gene_names <- str_match(TGF_beta_signaling_pathway, "\\(RefSeq\\)\\s*(.*?)\\s*,")
#remove NA's
gene_names <- na.omit(gene_names) 

# select column with gene names
gene_names <- as.data.frame(gene_names) #needs to be a df otherwise function select will not work
gene_names<-gene_names %>% select(V2)

#write the gene list to csv
write.csv(x=gene_names, file="TGF_beta_signaling_pathway_gene_names")



######### mevalonate pathway #########
### terpenoid_backbone_biosynthesis ###
#read file as data frame
terpenoid_backbone_biosynthesis <- read.delim2(file = "terpenoid_backbone_biosynthesis.txt", header = FALSE, sep = "\n",dec = ",")
# read as character vector
terpenoid_backbone_biosynthesis <- sapply(terpenoid_backbone_biosynthesis, as.character) 

#regex for gene names
gene_names <- str_match(terpenoid_backbone_biosynthesis, "\\(RefSeq\\)\\s*(.*?)\\s*,")
#remove NA's
gene_names <- na.omit(gene_names) 

# select column with gene names
gene_names <- as.data.frame(gene_names) #needs to be a df otherwise function select will not work
gene_names<-gene_names %>% select(V2)

#write the gene list to csv
write.csv(x=gene_names, file="terpenoid_backbone_biosynthesis_gene_names")


