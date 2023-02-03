### Ejercicio ###

## Descargamos datos
speaqeasy_data <- file.path(tempdir(), "rse_speaqeasy.RData")
download.file("https://github.com/LieberInstitute/SPEAQeasy-example/blob/master/rse_speaqeasy.RData?raw=true", speaqeasy_data, mode = "wb")

## Cargamos SummarizedExperiment
library("SummarizedExperiment")
# Cargando objeto rse_gene
load(speaqeasy_data, verbose = TRUE)

## Exploramos el objeto
rse_gene

## Eliminemos el diagnosis "Other" porque no tiene información
rse_gene$PrimaryDx <- droplevels(rse_gene$PrimaryDx)
table(rse_gene$PrimaryDx)


### ¿Hay diferencias en totalAssignedGene o mitoRate entre los grupos de diagnosis (PrimaryDx)? ###

## Exploramos PrimaryDx
table(rse_gene$PrimaryDx)

## Exploramos las diferencias entre grupos de diagnosis para varias variables

# totalAssignedGene
with(colData(rse_gene), tapply(totalAssignedGene, PrimaryDx, summary))

# mitoRate
with(colData(rse_gene), tapply(mitoRate, PrimaryDx, summary))


### Grafica la expresión de SNAP25 para cada grupo de diagnosis ###

# Explora el objeto rse de forma interactiva
library("iSEE")
iSEE::iSEE(rse_gene)

### Sugerir un modelo estadístico para usar en un análisis de expresión diferencial ###

