#' Data processing and summarization of peptide-level quantification to PTM and protein level quantification
#' 
#' Function to perform data processing and summarization on an experiment 
#' targeting post-translational modifications. Performs normalization, missing 
#' value imputation, feature selection, and summarization. Can optionally take 
#' an additional global protein quantification experiment for protein-level 
#' correction of PTM changes. Can take either label free or tandem mass tag 
#' (TMT) labeled data.
#' 
#' @export
#' @importFrom MSstats dataProcess
#' @importFrom MSstatsTMT proteinSummarization
#' 
#' @param data Name of the output of MSstatsPTM converter function or 
#' peptide-level quantified data from other tools. It should be a list 
#' containing one or two data tables, named PTM and PROTEIN for modified and 
#' unmodified datasets. The list must at least contain the PTM dataset, however
#' the PROTEIN dataset is optional.
#' @param ptm_label_type Indicator of labeling type for PTM dataset. Must be one of
#' `LF` or `TMT`
#' @param protein_label_type Indicator of labeling type for PROTEIN dataset. Must be 
#' one of `LF` or `TMT`
#' @param MBimpute_ptm Missing value imputation on the feature-level for the PTM
#' dataset. TRUE (default) imputes missing values by Accelated 
#' failure model. FALSE uses minimum value to impute the missing value for each 
#' peptide precursor ion.
#' @param MBimpute_protein Missing value imputation on the feature-level for the
#' PROTEIN dataset. TRUE (default) imputes missing values by Accelated 
#' failure model. FALSE uses minimum value to impute the missing value for each 
#' peptide precursor ion.
#' @param ... Additional parameters passed to either the 
#' \link[MSstats]{dataProcess} function from `MSstats` or the 
#' \link[MSstatsTMT]{proteinSummarization} function from `MSstatsTMT`.
dataProcessPTM = function(data, 
                          ptm_label_type,
                          protein_label_type,
                          MBimpute_ptm=FALSE,
                          MBimpute_protein=TRUE, 
                          use_log_file=TRUE, 
                          append=FALSE,
                          verbose=TRUE, 
                          log_file_path=NULL,
                          ...){
  
  ## Start log
  if (is.null(log_file_path) & use_log_file == TRUE){
    time_now = Sys.time()
    path = paste0("MSstatsPTM_log_", gsub("[ :\\-]", "_", time_now),
                  ".log")
    file.create(path)
  } else {
    path = log_file_path
  }
  
  ars = list(...)
  
  potential_lf_ars = c("logTrans", "normalization", "nameStandards", 
                       "featureSubset", "remove_uninformative_feature_outlier", 
                       "min_feature_count", "n_top_feature", "summaryMethod", 
                       "equalFeatureVar", "censoredInt", "remove50missing", 
                       "fix_missing", "maxQuantileforCensored", "numberOfCores")
  potential_tmt_ars = c("method", "global_norm", "reference_norm", 
                        "remove_norm_channel", "remove_empty_channel", 
                        "maxQuantileforCensored")
  
  lf_ars = intersect(names(ars), potential_lf_ars)
  tmt_ars = intersect(names(ars), potential_tmt_ars)
  
  ## Determine if PTM should be adjusted for protein level.
  if (!is.null(data$PROTEIN)){
    adj_protein = TRUE
  } else{
    adj_protein = FALSE
  }
  
  if (ptm_label_type == "LF"){
    
    MSstatsLogsSettings(use_log_file, append,
                        verbose, log_file_path = path)
    
    summarized_ptm = do.call(
      dataProcess, c(list(raw=data$PTM), 
                     ars[lf_ars],
                     list(MBimpute=MBimpute_ptm))
    )
  } else if (ptm_label_type == "TMT"){
    
    MSstatsLogsSettings(use_log_file, append,
                        verbose, log_file_path = path,
                        pkg_name = "MSstatsTMT")
    
    summarized_ptm = do.call(
      proteinSummarization, c(list(data=data$PTM), 
                              ars[potential_tmt_ars],
                              list(MBimpute=MBimpute_ptm))
    )

  }
  if (adj_protein){
    if (protein_label_type == "LF"){
      
      MSstatsLogsSettings(use_log_file, append,
                          verbose, log_file_path = path)
      
      summarized_protein = do.call(
        dataProcess, c(list(raw=data$PROTEIN), 
                       ars[lf_ars],
                       list(MBimpute=MBimpute_protein))
      )
    } else if (protein_label_type == "TMT"){
      
      MSstatsLogsSettings(use_log_file, append,
                          verbose, log_file_path = path,
                          pkg_name = "MSstatsTMT")
      
      summarized_protein = do.call(
        proteinSummarization, c(list(data=data$PROTEIN), 
                                ars[potential_tmt_ars],
                                list(MBimpute=MBimpute_protein))
      )
    }
  }
  
  return(
    list(
      PTM=summarized_ptm,
      PROTEIN=summarized_protein
      )
    )
  
}