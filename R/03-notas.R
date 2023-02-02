### Modelos estadísticos en R ###

## model.matrix
mat <- with(trees, model.matrix(log(Volume) ~ log(Height) + log(Girth)))
mat

colnames(mat)

## Summary
summary(lm(log(Volume) ~ log(Height) + log(Girth), data = trees))

### ExploreModelMatrix ###
## Es un paquete de Bioconductor que nos ayuda a entender los modelos estadísticos
## que estamos usando gracias a visualizaciones

## Datos de ejemplo
(sampleData <- data.frame(
  genotype = rep(c("A", "B"), each = 4),
  treatment = rep(c("ctrl", "trt"), 4)
))

## Creemos las imágenes usando ExploreModelMatrix
vd <- ExploreModelMatrix::VisualizeDesign(
  sampleData = sampleData,
  designFormula = ~ genotype + treatment,
  textSizeFitted = 4
)

## Veamos las imágenes con cowplot
cowplot::plot_grid(plotlist = vd$plotlist)

## Usaremos shiny otra ves
app <- ExploreModelMatrix(
  sampleData = sampleData,
  designFormula = ~ genotype + treatment
)
if (interactive()) shiny::runApp(app)


### Ejemplo 2 ###

(sampleData <- data.frame(
  Response = rep(c("Resistant", "Sensitive"), c(12, 18)),
  Patient = factor(rep(c(1:6, 8, 11:18), each = 2)),
  Treatment = factor(rep(c("pre","post"), 15)),
  ind.n = factor(rep(c(1:6, 2, 5:12), each = 2))))

vd <- VisualizeDesign(
  sampleData = sampleData,
  designFormula = ~ Response + Response:ind.n + Response:Treatment,
  textSizeFitted = 3
)

cowplot::plot_grid(plotlist = vd$plotlist, ncol = 1)

### Ejemplo 3 ###

(sampleData = data.frame(
  condition = factor(rep(c("ctrl_minus", "ctrl_plus",
                           "ko_minus", "ko_plus"), 3)),
  batch = factor(rep(1:6, each = 2))))

# 0 + se usa para no tener intercepto (no quieres algo específico como referencia)
vd <- VisualizeDesign(sampleData = sampleData,
                      designFormula = ~ 0 + batch + condition,
                      textSizeFitted = 4, lineWidthFitted = 20,
                      dropCols = "conditionko_minus")

cowplot::plot_grid(plotlist = vd$plotlist, ncol = 1)


###  Datos de SRP045638 con recount 3 ###

library("recount3")

options(recount3_url = "https://recount-opendata.s3.amazonaws.com/recount3/release")

human_projects <- available_projects()

rse_gene_SRP045638 <- create_rse(
  subset(
    human_projects,
    project == "SRP045638" & project_type == "data_sources"
  )
)

assay(rse_gene_SRP045638, "counts") <- compute_read_counts(rse_gene_SRP045638)

# Una vez descargados y con los números de lecturas podemos usar expand_sra_attributes().
# Sin embargo, tenemos un problema con estos datos.

rse_gene_SRP045638$sra.sample_attributes[1:3]

# Vamos a intentar resolverlo eliminando información que está presente solo en ciertas muestras.

rse_gene_SRP045638$sra.sample_attributes <- gsub("dev_stage;;Fetal\\|", "", rse_gene_SRP045638$sra.sample_attributes)
rse_gene_SRP045638$sra.sample_attributes[1:3]

# Repetimos el código de ayer
rse_gene_SRP045638 <- expand_sra_attributes(rse_gene_SRP045638)

colData(rse_gene_SRP045638)[
  ,
  grepl("^sra_attribute", colnames(colData(rse_gene_SRP045638)))
]


## Pasar de character a nuemric o factor
rse_gene_SRP045638$sra_attribute.age <- as.numeric(rse_gene_SRP045638$sra_attribute.age)
rse_gene_SRP045638$sra_attribute.disease <- factor(rse_gene_SRP045638$sra_attribute.disease)
rse_gene_SRP045638$sra_attribute.RIN <- as.numeric(rse_gene_SRP045638$sra_attribute.RIN)
rse_gene_SRP045638$sra_attribute.sex <- factor(rse_gene_SRP045638$sra_attribute.sex)

## Resumen de las variables de interés
summary(as.data.frame(colData(rse_gene_SRP045638)[
  ,
  grepl("^sra_attribute.[age|disease|RIN|sex]", colnames(colData(rse_gene_SRP045638)))
]))

### Encontraremos diferencias entre muestra prenatalas vs postnatales ###

rse_gene_SRP045638$prenatal <- factor(ifelse(rse_gene_SRP045638$sra_attribute.age < 0, "prenatal", "postnatal"))
table(rse_gene_SRP045638$prenatal)

## http://rna.recount.bio/docs/quality-check-fields.html
rse_gene_SRP045638$assigned_gene_prop <- rse_gene_SRP045638$recount_qc.gene_fc_count_all.assigned / rse_gene_SRP045638$recount_qc.gene_fc_count_all.total
summary(rse_gene_SRP045638$assigned_gene_prop)

### Graficamos ###
with(colData(rse_gene_SRP045638), plot(assigned_gene_prop, sra_attribute.RIN))

# Diferencias entre los grupos
with(colData(rse_gene_SRP045638), tapply(assigned_gene_prop, prenatal, summary))

## Guardemos nuestro objeto entero por si luego cambiamos de opinión
rse_gene_SRP045638_unfiltered <- rse_gene_SRP045638

## Eliminemos a muestras malas
hist(rse_gene_SRP045638$assigned_gene_prop)

table(rse_gene_SRP045638$assigned_gene_prop < 0.3)

# Nos quedamos con los datos buenos
rse_gene_SRP045638 <- rse_gene_SRP045638[, rse_gene_SRP045638$assigned_gene_prop > 0.3]


## Calculemos los niveles medios de expresión de los genes en nuestras
## muestras.
## Ojo: en un análisis real probablemente haríamos esto con los RPKMs o CPMs
## en vez de las cuentas.
gene_means <- rowMeans(assay(rse_gene_SRP045638, "counts"))
summary(gene_means)

## Eliminamos genes
rse_gene_SRP045638 <- rse_gene_SRP045638[gene_means > 0.1, ]

## Dimensiones finales
dim(rse_gene_SRP045638)

## Porcentaje de genes que retuvimos
round(nrow(rse_gene_SRP045638) / nrow(rse_gene_SRP045638_unfiltered) * 100, 2)

### Normalización  ###

library("edgeR") # BiocManager::install("edgeR", update = FALSE)
dge <- DGEList(
  counts = assay(rse_gene_SRP045638, "counts"),
  genes = rowData(rse_gene_SRP045638)
)
dge <- calcNormFactors(dge)


### Exploramos la relación entre nuestras variables: Expresión diferencial ###

library("ggplot2")

ggplot(as.data.frame(colData(rse_gene_SRP045638)), aes(y = assigned_gene_prop, x = prenatal)) +
  geom_boxplot() +
  theme_bw(base_size = 20) +
  ylab("Assigned Gene Prop") +
  xlab("Age Group")

