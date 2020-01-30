###############################################################################
####                             OMIM DATA                                 ####
###############################################################################

## Cargamos packages
library(biomaRt)

## Cargamos los datos de OMIM
omim <- read.table("omim_data_20200123.txt", header = FALSE, sep = "\t")

## Cargamos el nombre de los genes que queremos obtener los datos del OMIM
genes <- read.table("genes.txt", header = FALSE)
genes <- as.character(unlist(genes))

## Creamos un nuevo objeto con los ID de Ensembl de los genes que queremos
omim_genes_f <- omim$V4 %in% genes
ens_id <- omim[omim_genes_f,5]
omim_genes <- omim[omim_genes_f,]
rownames(omim_genes) <- omim_genes$V4
## Obtenemos los datos de Ensembl con el package BiomaRt
biomart_uses <- useMart(biomart="ensembl", dataset="hsapiens_gene_ensembl")

OMIM_data_genes <- getBM(attributes=c("hgnc_symbol","mim_morbid_description", 
                                      "mim_morbid_accession"),
                          filter="ensembl_gene_id",
                          values=ens_id,
                          mart=biomart_uses,
                          uniqueRows=TRUE)

location_ens <- getBM(attributes=c("chromosome_name","band"),
                      filter="ensembl_gene_id",
                      values=ens_id,
                      mart=biomart_uses,
                      uniqueRows=TRUE)

OMIM_data_genes$location <- apply(location_ens, 1, function(x){
  paste(x[1],x[2], sep = "")
})

rownames(OMIM_data_genes) <- OMIM_data_genes[,1]

## De la columna de mim morbid description solo nos quedamos con lo que nos 
# interesa (lo que hay antes del primer ";")
description_list <- lapply(OMIM_data_genes$mim_morbid_description, function(x){
  y <- unlist(strsplit(x, ";"))
  out <- y[1]
  return(out)
})

OMIM_data_genes$mim_morbid_description <- do.call("rbind", description_list)

## Unimos toda la información que nos interesa en un data.frame para poder
# exportarlo en formato txt y poder abrir el archivo con excel.

omim_union <- cbind(OMIM_data_genes$hgnc_symbol, omim_genes$V1,
                    OMIM_data_genes$location,
                    OMIM_data_genes$mim_morbid_description,
                    OMIM_data_genes$mim_morbid_accession)
colnames(omim_union) <- c("Nom_gen_HGNC",	"Numero_gen_MIM", 
                          "Localitzacio_gen",
                          "Fenotip_asociat", "MIM_Fenotip") 

omim_union[,3] <- tolower(omim_union[,3])

## Guardamos el objeto omim_union en un archivo txt 
write.table(omim_union, file = "dataomim_gene.txt", sep = "\t", quote=FALSE, row.names = FALSE)

