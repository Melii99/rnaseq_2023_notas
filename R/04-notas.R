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


# También es posible hacer graficas nosotros mismos. Aquí les muestro una posible respuesta
# con ggplot2
library("ggplot2")
ggplot(
  as.data.frame(colData(rse_gene)),
  aes(y = totalAssignedGene, group = PrimaryDx, x = PrimaryDx)
) +
  geom_boxplot() +
  theme_bw(base_size = 20) +
  xlab("Diagnosis")

### Grafica la expresión de SNAP25 para cada grupo de diagnosis ###

# Explora el objeto rse de forma interactiva
library("iSEE")
iSEE::iSEE(rse_gene)

# Respuesta alternativa con ggplot2

## Encontremos el gene SNAP25
rowRanges(rse_gene)
## En este objeto los nombres de los genes vienen en la variable "Symbol"
i <- which(rowRanges(rse_gene)$Symbol == "SNAP25")
i
## Para graficar con ggplot2, hagamos un pequeño data.frame
df <- data.frame(
  expression = assay(rse_gene)[i, ],
  Dx = rse_gene$PrimaryDx
)
## Ya teniendo el pequeño data.frame, podemos hacer la gráfica
ggplot(df, aes(y = log2(expression + 0.5), group = Dx, x = Dx)) +
  geom_boxplot() +
  theme_bw(base_size = 20) +
  xlab("Diagnosis") +
  ylab("SNAP25: log2(x + 0.5)")


### Sugerir un modelo estadístico para usar en un análisis de expresión diferencial ###

