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
with(sampleData, model.matrix(~ genorype + treatment))

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
