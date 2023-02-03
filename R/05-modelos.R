## Model matrix
#Matriz que define el modelo estadístico.

mat <- with(trees, model.matrix(log(Volume) ~ log(Height) + log(Girth)))
mat

#obtener  de una regresion lineal.
summary(lm(log(Volume) ~ log(Height) + log(Girth), data = trees))
#primer numero es el intercepto de la variable y cuando x == 0
## output
##  Coefficients:
##              Estimate Std. Error t value Pr(>|t|)
##  (Intercept) -6.63162    0.79979  -8.292 5.06e-09
##                este

################ ExploreModelMatrix
## Datos de ejemplo
(sampleData <- data.frame(
  genotype = rep(c("A", "B"), each = 4),
  treatment = rep(c("ctrl", "trt"), 4)
))

# Para verlo con unos y ceros
with(sampleData, model.matrix(~ genotype + treatment))

## Creemos las imágenes usando ExploreModelMatrix
vd <- ExploreModelMatrix::VisualizeDesign(
  sampleData = sampleData,
  designFormula = ~ genotype + treatment,
  textSizeFitted = 4
)

## Veamos las imágenes
cowplot::plot_grid(plotlist = vd$plotlist)

#para isntalar
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("ExploreModelMatrix")

## Usaremos shiny otra ves
app <- ExploreModelMatrix(
  sampleData = sampleData,
  designFormula = ~ genotype + treatment
)
if (interactive()) shiny::runApp(app)

#Interpreta ResponseResistant.Treatmentpre del ejercicio 2.
#Puede ser útil tomar un screenshot (captura de pantalla)
#y anotarla con líneas de colores.
#Si haces eso, puedes incluir la imagen en tus notas.

###Respuesta:
#Tenemos que restar pre menos post tratamiento.
# ResponseResistant:Treatmentpre = pre - post
# Si tenemos -7 dirías que hay más expresión después del tratamiento.

# ¿Por qué es clave el 0 al inicio de la fórmula en el ejercicio 3?
# Para que no tome en cuenta el intercepto y puedes saber
#cuales son los valores delos batch solos.

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

#Vemos que hay datos que no estan alineados
rse_gene_SRP045638$sra.sample_attributes[1:3]

#Editar los datos para que esten alineados todas las tablas.
rse_gene_SRP045638$sra.sample_attributes <- gsub("dev_stage;;Fetal\\|", "", rse_gene_SRP045638$sra.sample_attributes)
rse_gene_SRP045638$sra.sample_attributes[1:3]

#Ver los datos
rse_gene_SRP045638 <- expand_sra_attributes(rse_gene_SRP045638)

colData(rse_gene_SRP045638)[
  ,
  grepl("^sra_attribute", colnames(colData(rse_gene_SRP045638)))
]

## Pasar de character a numeric o factor
rse_gene_SRP045638$sra_attribute.age <- as.numeric(rse_gene_SRP045638$sra_attribute.age)
rse_gene_SRP045638$sra_attribute.disease <- factor(tolower(rse_gene_SRP045638$sra_attribute.disease))
rse_gene_SRP045638$sra_attribute.RIN <- as.numeric(rse_gene_SRP045638$sra_attribute.RIN)
rse_gene_SRP045638$sra_attribute.sex <- factor(rse_gene_SRP045638$sra_attribute.sex)

#
rse_gene_SRP045638$prenatal <- factor(ifelse(rse_gene_SRP045638$sra_attribute.age < 0, "prenatal", "postnatal"))

#
rse_gene_SRP045638$prenatal <- factor(ifelse(rse_gene_SRP045638$sra_attribute.age < 0, "prenatal", "postnatal"))
table(rse_gene_SRP045638$prenatal)

#proporcion de muestras que se alinearon con genes
rse_gene_SRP045638$assigned_gene_prop <- rse_gene_SRP045638$recount_qc.gene_fc_count_all.assigned / rse_gene_SRP045638$recount_qc.gene_fc_count_all.total
summary(rse_gene_SRP045638$assigned_gene_prop)

## RIN : RNA integrity number. mientras mas alto menos rnas rotos.
#Correlacion de RIN con la proporcion de lecturas asignadas a genes
with(colData(rse_gene_SRP045638), plot(assigned_gene_prop, sra_attribute.RIN))

rse_gene_SRP045638_unfiltered <- rse_gene_SRP045638

## Eliminemos a muestras malas
hist(rse_gene_SRP045638$assigned_gene_prop)

rse_gene_SRP045638 <- rse_gene_SRP045638[, rse_gene_SRP045638$assigned_gene_prop > 0.3]

## Calculemos los niveles medios de expresión de los genes en nuestras
## muestras.
## Ojo: en un análisis real probablemente haríamos esto con los
# RPKMs (reads) o CPMs ()
## en vez de las cuentas.
gene_means <- rowMeans(assay(rse_gene_SRP045638, "counts"))
summary(gene_means)

## Eliminamos genes con casi nada de expresion
rse_gene_SRP045638 <- rse_gene_SRP045638[gene_means > 0.1, ]

## Dimensiones finales
dim(rse_gene_SRP045638)

library("edgeR") # BiocManager::install("edgeR", update = FALSE)
dge <- DGEList(
  counts = assay(rse_gene_SRP045638, "counts"),
  genes = rowData(rse_gene_SRP045638)
)
dge <- calcNormFactors(dge)

#### Expresion diferencial
#Definir el modelo estadistico.
library("ggplot2")
ggplot(as.data.frame(colData(rse_gene_SRP045638)), aes(y = assigned_gene_prop, x = prenatal)) +
  geom_boxplot() +
  theme_bw(base_size = 20) +
  ylab("Assigned Gene Prop") +
  xlab("Age Group")

mod <- model.matrix(~ prenatal + sra_attribute.RIN + sra_attribute.sex + assigned_gene_prop,
                    data = colData(rse_gene_SRP045638)
)
colnames(mod)

#
library("limma")
vGene <- voom(dge, mod, plot = TRUE)

#
eb_results <- eBayes(lmFit(vGene))

de_results <- topTable(
  eb_results,
  coef = 2,
  number = nrow(rse_gene_SRP045638),
  sort.by = "none"
)

#
head(de_results)
