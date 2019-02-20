#!/usr/local/bin/Rscript

if (!require("pacman")) install.packages("pacman")
library(pacman)
p_load(optparse)


option_list = list(
  make_option(c("-m", "--model_pred"), type="character", default=NULL, 
              help="model prediction CSV filename", metavar="file"),
  make_option(c("-t", "--model_type"), type="character", default=NULL, 
              help="Must be either RF or GB", metavar="character")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$model_pred)){
  print_help(opt_parser)
  stop("There is no model prediction file!.\n", call.=FALSE)
}

if (is.null(opt$model_type)){
    print_help(opt_parser)
    step("You must select a model type: GB (default) or RF.\n", call.=FALSE)
}

getScriptPath <- function(){
    cmd.args <- commandArgs()
    m <- regexpr("(?<=^--file=).+", cmd.args, perl=TRUE)
    script.dir <- dirname(regmatches(cmd.args, m))
    if(length(script.dir) == 0) stop("can't determine script dir: please call the script with Rscript")
    if(length(script.dir) > 1) stop("can't determine script dir: more than one '--file' argument detected")
    return(script.dir)
}

script_dir <- getScriptPath()
#cat("file path relative to current directory", script_dir, "\n")
script_dir <- substr(script_dir,1,nchar(script_dir)-4)
#cat("path to crema relative to current directory", script_dir, "\n")
file_path_file <- file.path(script_dir, "log_reg_models", "logit_models.RData")
#cat("path to logit_models.RData", file_path_file, "\n")
load(file_path_file)

predict_new <- function(model_output_csv, gene_names, model_name, model_type){
    model <- read.csv(model_output_csv, row.names=1)
    geneNames <- rownames(model)
    geneNames <- make.unique(geneNames)
    #geneNames <- make.unique(as.character(geneNames[,1]))
    #modelPred <- model[,c(1,3,5,7,9,11,13,15)]
    colnames(model) <- paste0(model_type, c(1:8))
    #rownames(modelPred) <- geneNames[,1]
    modelPred <- as.data.frame(apply(model, 2, as.factor))
    rownames(modelPred) <- geneNames
    model_logitPred <- predict(model_name, newdata=modelPred, type="response")
    return(model_logitPred)
}


if (opt$model_type == "RF"){
	model_name <- RF_logit
}

if (opt$model_type == "GB"){
	model_name <- GB_logit
}


new_predictions <- predict_new(opt$model_pred, opt$gene_name, model_name, opt$model_type)
file_name <- "ensemble_logreg_pred.csv"
write.csv(new_predictions, file=file_name)


