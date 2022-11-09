library(AnnotationHub)
library(ensembldb)
library(AnnotationForge)
# start an AnnotationHub instance/connection.
ah <- AnnotationHub()
# query for availabel Rat Ensembl databases
EnsDb.rat <- query(ah, c("EnsDb", "Rattus norvegicus"))
EnsDb.rat
edb <- ah[["AH75088"]]
setwd("/work/LAS/geetu-lab/hhvu/project3_scATAC/scATAC-seq-analysis/")
#copy databse from the cache to working dir
file.copy(AnnotationHub::cache(ah["AH75088"]), "./EnsDb.rat.sqlite")
# now make it a package. Change name and email accordingly
makeEnsembldbPackage("EnsDb.rat.sqlite", version="0.0.1",
                     maintainer = "First Name <first.last@mydomain.com>",
                     author = "First Name <first.last@mydomain.com>",
                     destDir=".", license="Artistic-2.0")

# install package in R.  
# note modifed name of created directory (not EnsDb.rat... but EnsDb.Rnorvegicus.v100)
install.packages("./EnsDb.Rnorvegicus.v98", type = "source", repos = NULL)

#check: 
library(EnsDb.Rnorvegicus.v98)
EnsDb.Rnorvegicus.v98

columns(EnsDb.Rnorvegicus.v98)
#[1] "DESCRIPTION"         "ENTREZID"            "EXONID"             
#<<snip>>      
#[37] "UNIPROTID"           "UNIPROTMAPPINGTYPE" 