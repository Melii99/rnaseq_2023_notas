## Lets build our first SummarizedExperiment object
library("SummarizedExperiment")
## ?SummarizedExperiment

## De los ejemplos en la ayuda oficial

## Creamos los datos para nuestro objeto de tipo SummarizedExperiment
## para 200 genes a lo largo de 6 muestras
nrows <- 200
ncols <- 6
## Números al azar de cuentas
set.seed(20210223)
counts <- matrix(runif(nrows * ncols, 1, 1e4), nrows)
## Información de nuestros genes
rowRanges <- GRanges(
  rep(c("chr1", "chr2"), c(50, 150)),
  IRanges(floor(runif(200, 1e5, 1e6)), width = 100),
  strand = sample(c("+", "-"), 200, TRUE),
  feature_id = sprintf("ID%03d", 1:200)
)
names(rowRanges) <- paste0("gene_", seq_len(length(rowRanges)))
## Información de nuestras muestras
colData <- DataFrame(
  Treatment = rep(c("ChIP", "Input"), 3),
  row.names = LETTERS[1:6]
)
## Juntamos ahora toda la información en un solo objeto de R
rse <- SummarizedExperiment(
  assays = SimpleList(counts = counts),
  rowRanges = rowRanges,
  colData = colData
)

## Exploremos el objeto resultante
rse

### Explorando el objeto ###

## Número de genes y muestras
dim(rse)

## IDs de nuestros genes y muestras
dimnames(rse)

## Nombres de tablas de cuentas que tenemos (RPKM, CPM, counts, logcounts, etc)
assayNames(rse)

## El inicio de nuestra tabla de cuentas
head(assay(rse))

## Información de los genes en un objeto de Bioconductor
rowRanges(rse)

## Tabla con información de los genes
rowData(rse) # es idéntico a 'mcols(rowRanges(rse))'

## Tabla con información de las muestras
colData(rse)

### Ejercicio ###

## Comando 1
rse[1:2, ] # Se imprime el primer y segundo renglón del objeto y todas las columnas

## Comando 2
rse[, c("A", "D", "F")] # Se imprimen todos los renglones del objeto y las columnas A,D y F

### iSEE ###

## Explora el objeto rse de forma interactiva
library("iSEE")
iSEE::iSEE(rse)

### Ejercicio con spatialLIBD ###

## Descarguemos unos datos de spatialLIBD
sce_layer <- spatialLIBD::fetch_data("sce_layer")

## Revilsemos el tamaño de este objeto
obstr::obj_size(sce_layer)

sce_layer
rowData(sce_layer)
colData(sce_layer)

# Al igual que nuestro objeto rse podemos usar iSEE::iSEE() para explorar los datos.
iSEE::iSEE(sce_layer)


### recount3 ###

## Load recount3 R package
library("recount3")

## Cambiaremos el URL de recount3 a Amazon (AWS)
## para evitar problems de TLS versión 1.2 en
## Windows (luego tiene una versión más vieja).

## Primero vemos el URL actual
getOption(
  "recount3_url",
  "http://duffel.rail.bio/recount3"
)

## Ahora si lo cambiamos
options(recount3_url = "https://recount-opendata.s3.amazonaws.com/recount3/release")

## Y vemos que ya funcionó el cambio.
getOption(
  "recount3_url",
  "http://duffel.rail.bio/recount3"
)

## Revisemos todos los proyectos con datos de humano en recount3
human_projects <- available_projects()

## Encuentra tu proyecto de interés. Aquí usaremos
## SRP009615 de ejemplo
proj_info <- subset(
  human_projects,
  project == "SRP009615" & project_type == "data_sources"
)
## Crea un objeto de tipo RangedSummarizedExperiment (RSE)
## con la información a nivel de genes
rse_gene_SRP009615 <- create_rse(proj_info)

## Explora el objeto RSE
rse_gene_SRP009615

## Explora los proyectos disponibles de forma interactiva
proj_info_interactive <- interactiveDisplayBase::display(human_projects)
## Selecciona un solo renglón en la tabla y da click en "send".

## Aquí verificamos que solo seleccionaste un solo renglón.
stopifnot(nrow(proj_info_interactive) == 1)
## Crea el objeto RSE
rse_gene_interactive <- create_rse(proj_info_interactive)

## Convirtamos las cuentas por nucleotido a cuentas por lectura
## usando compute_read_counts().
## Para otras transformaciones como RPKM y TPM, revisa transform_counts().
assay(rse_gene_SRP009615, "counts") <- compute_read_counts(rse_gene_SRP009615)

## Para este estudio en específico, hagamos más fácil de usar la
## información del experimento
rse_gene_SRP009615 <- expand_sra_attributes(rse_gene_SRP009615)
colData(rse_gene_SRP009615)[
  ,
  grepl("^sra_attribute", colnames(colData(rse_gene_SRP009615)))
]

### Ejercicio iSEE imágen ###
iSEE::iSEE(rse_gene_SRP009615)
