library("recount3")
## Revisemos todos los proyectos con datos de humano en recount3
available_projects(recount3_url = "https://recount-opendata.s3.amazonaws.com/recount3/release")
options(recount3_url = "https://recount-opendata.s3.amazonaws.com/recount3/release")

human_projects <- available_projects()
## Encuentra tu proyecto de interés. Aquí usaremos
## SRP009615 de ejemplo
proj_info <- subset(
  human_projects,
  project == "SRP009615" & project_type == "data_sources"
)
## Crea un objetio de tipo RangedSummarizedExperiment (RSE)
## con la información a nivel de genes
rse_gene_SRP009615 <- create_rse(proj_info)
## Explora el objeto RSE
rse_gene_SRP009615


## ----"interactive_display", eval = FALSE------------
## ## Explora los proyectos disponibles de forma interactiva
## proj_info_interactive <- interactiveDisplayBase::display(human_projects)
## ## Selecciona un solo renglón en la tabla y da click en "send".
##
## ## Aquí verificamos que solo seleccionaste un solo renglón.
## stopifnot(nrow(proj_info_interactive) == 1)
## ## Crea el objeto RSE
## rse_gene_interactive <- create_rse(proj_info_interactive)

# base pair coverage es el numero de alineamientos que coincidieron
# en esa base. Las lecturas que cubrieron esa base.


## ----"tranform_counts"------------------------------
## Convirtamos las cuentas por nucleotido a cuentas por lectura
## usando compute_read_counts().
## Para otras transformaciones como RPKM y TPM, revisa transform_counts().
assay(rse_gene_SRP009615, "counts") <- compute_read_counts(rse_gene_SRP009615)


## ----"expand_attributes"----------------------------
## Para este estudio en específico, hagamos más fácil de usar la
## información del experimento
rse_gene_SRP009615 <- expand_sra_attributes(rse_gene_SRP009615)
colData(rse_gene_SRP009615)[
  ,
  grepl("^sra_attribute", colnames(colData(rse_gene_SRP009615)))
]
