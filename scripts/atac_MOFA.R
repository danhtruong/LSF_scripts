#initialize libraries

library(reticulate)
py_discover_config(required_module = 'mofapy2')
py_module_available('mofapy2')
library(Seurat)
library(MOFA2)
library(cisTopic)

#load cistopic data
cisTopicObject <- readRDS('2020-06-03 cistopicobject.rds')

cistopics_embeddings <- modelMatSelection(
  cisTopicObject, 
  target = "cell", 
  method = "Z-score"
)

#load chromvar data
motif_assay <- readRDS('2020-06-07 MTL_ATAC_chromvar_assay.rds')

#Load metadata
MTL_ATAC.subset <- readRDS('2020-06-07 MTL_ATAC.metadata.rds')

#create MOFA2
mofa <- create_mofa(list(
  "ATAC" = cistopics_embeddings,
  "motifs" = motif_assay@data),
  groups = MTL_ATAC.subset$Condition)

# Model options
model_opts <- get_default_model_options(mofa)
model_opts$num_factors <- 10

# Data options
data_opts <- get_default_data_options(mofa)
data_opts$scale_views

# Training options
train_opts <- get_default_training_options(mofa)
train_opts$seed <- 42

# Prepare MOFA object
mofa <- prepare_mofa(
  object = mofa,
  model_options = model_opts,
  training_options = train_opts
)

#run the model
outfile = paste0(getwd(),"/atac_cistopic_chromvar_model.hdf5")
#model <- load_model(outfile)
model <- run_mofa(mofa, outfile = outfile)