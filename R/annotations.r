# ===========================================================================
# R/annotations.R — User-customizable annotation helpers and loader
# ===========================================================================

#' Load Built-in Curated or Freshly Scored Biological Annotations
#'
#' Retrieves either the frozen curated benchmark annotation matrix (Supplementary Table S1.21)
#' or allows building fresh scores using automated deterministic rules from raw database queries.
#'
#' @param key_type Character, identifier column to use: "entrez" (ID), "ensembl" (ENSEMBL_GENE_ID), or "symbol" (GeneName).
#' @param version Character, either \code{"curated"} (default, the validated benchmark matrix used in the manuscript)
#'   or \code{"raw_scored"} (when \code{raw_df} is supplied, recalculates using deterministic scoring rules).
#' @param raw_df Optional data frame containing raw DAVID query outputs. Only used if \code{version = "raw_scored"}.
#' @param col_mapping Optional named list mapping column names in \code{raw_df} when \code{version = "raw_scored"}.
#' @return A data frame containing gene identifiers and biological annotation columns.
#' @export
#' @examples
#' # Option 1: Load the benchmark curated frozen matrix (Default)
#' annot_curated <- load_builtin_annotations("entrez", version = "curated")
#' head(annot_curated[, 1:5])
#'
#' # Option 2: Score raw database query table using reproducible rules
#' raw_toy <- data.frame(
#'   ID = c("1", "2"),
#'   tissue = c("Fetal brain, cortex", "Liver"),
#'   disease = c("Schizophrenia", "None"),
#'   pathway = c("Dopaminergic synapse", NA),
#'   tfbs_count = c(3, 0),
#'   process = c("Synaptic transmission", NA),
#'   interaction = c("Dopamine receptor binding", NA)
#' )
#' annot_fresh <- load_builtin_annotations("entrez", version = "raw_scored", raw_df = raw_toy)
#' head(annot_fresh)
load_builtin_annotations <- function(key_type = c("entrez", "ensembl", "symbol"),
                                     version = c("curated", "raw_scored"),
                                     raw_df = NULL,
                                     col_mapping = list(
                                       tissue = "tissue",
                                       disease = "disease",
                                       pathway = "pathway",
                                       tfbs = "tfbs_count",
                                       process = "process",
                                       interaction = "interaction"
                                     )) {
  key_type <- match.arg(key_type)
  version  <- match.arg(version)

  if (version == "raw_scored") {
    if (is.null(raw_df)) {
      stop("When version = 'raw_scored', 'raw_df' containing raw database queries must be provided.", call. = FALSE)
    }
    id_name <- if (key_type == "entrez") "ID" else if (key_type == "ensembl") "ENSEMBL_GENE_ID" else "GeneName"
    if (!id_name %in% colnames(raw_df)) {
      stop(sprintf("Identifier column '%s' required for key_type = '%s' in raw_df.", id_name, key_type), call. = FALSE)
    }
    return(score_david_annotations(raw_df, id_col = id_name, col_mapping = col_mapping))
  }

  utils::data("psychiatric_annotations", package = "fdrregGenomics", envir = environment())
  df <- get("psychiatric_annotations", envir = environment())
  
  if (key_type == "entrez") {
    return(df[, !colnames(df) %in% c("ENSEMBL_GENE_ID", "GeneName"), drop = FALSE])
  } else if (key_type == "ensembl") {
    cols <- setdiff(colnames(df), c("ID", "GeneName"))
    return(df[, c("ENSEMBL_GENE_ID", setdiff(cols, "ENSEMBL_GENE_ID")), drop = FALSE])
  } else {
    cols <- setdiff(colnames(df), c("ID", "ENSEMBL_GENE_ID"))
    return(df[, c("GeneName", setdiff(cols, "GeneName")), drop = FALSE])
  }
}

#' Prepare and Combine Custom User Annotations
#'
#' Provides a flexible helper to merge, format, and validate user-supplied
#' annotation data (e.g., custom pathway scores, single-cell markers, epigenetic peaks)
#' with optional baseline annotations.
#'
#' @param user_annotations A data frame containing custom annotations.
#' @param id_col Character, name of the identifier column in \code{user_annotations}. Default is "ID".
#' @param base_annotations Optional data frame of existing/baseline annotations to merge with.
#' @param base_id_col Character, name of the identifier column in \code{base_annotations}. Defaults to \code{id_col}.
#' @param fill_na Numeric value to replace missing values (NA/Inf) with. Default is 0.
#' @return A processed data frame with the ID column and numeric feature columns, ready for FDRreg.
#' @export
#' @examples
#' custom_df <- data.frame(
#'   ID = c(1, 2, 9, 10),
#'   pathway_score = c(1.2, 0.5, 3.1, 0.0),
#'   is_hub = c(1, 0, 1, 0)
#' )
#' prep_annot <- prepare_custom_annotations(custom_df, id_col = "ID")
#' head(prep_annot)
prepare_custom_annotations <- function(user_annotations, id_col = "ID",
                                       base_annotations = NULL, base_id_col = NULL,
                                       fill_na = 0) {
  if (!is.data.frame(user_annotations)) {
    stop("'user_annotations' must be a data frame.", call. = FALSE)
  }
  if (!id_col %in% colnames(user_annotations)) {
    stop(sprintf("ID column '%s' not found in 'user_annotations'.", id_col), call. = FALSE)
  }

  user_annotations[[id_col]] <- as.character(user_annotations[[id_col]])

  if (!is.null(base_annotations)) {
    if (!is.data.frame(base_annotations)) {
      stop("'base_annotations' must be a data frame.", call. = FALSE)
    }
    if (is.null(base_id_col)) base_id_col <- id_col
    if (!base_id_col %in% colnames(base_annotations)) {
      stop(sprintf("Base ID column '%s' not found in 'base_annotations'.", base_id_col), call. = FALSE)
    }
    base_annotations[[base_id_col]] <- as.character(base_annotations[[base_id_col]])
    
    merged <- merge(base_annotations, user_annotations,
                    by.x = base_id_col, by.y = id_col, all = TRUE)
    colnames(merged)[colnames(merged) == base_id_col] <- id_col
    res <- merged
  } else {
    res <- user_annotations
  }

  id_vals <- res[[id_col]]
  feature_cols <- setdiff(colnames(res), id_col)
  
  if (length(feature_cols) == 0) {
    stop("No feature columns found in annotation data.", call. = FALSE)
  }

  num_mat <- as.matrix(res[, feature_cols, drop = FALSE])
  mode(num_mat) <- "numeric"
  
  num_mat[!is.finite(num_mat)] <- fill_na

  col_sds <- apply(num_mat, 2, function(x) stats::sd(x, na.rm = TRUE))
  bad_cols <- which(col_sds == 0 | is.na(col_sds))
  if (length(bad_cols) > 0) {
    warning(sprintf("[prepare_custom_annotations] Dropping %d zero-variance column(s): %s",
                    length(bad_cols), paste(colnames(num_mat)[bad_cols], collapse = ", ")),
            call. = FALSE)
    num_mat <- num_mat[, -bad_cols, drop = FALSE]
  }

  out <- data.frame(num_mat, check.names = FALSE)
  out[[id_col]] <- id_vals
  out <- out[, c(id_col, setdiff(colnames(out), id_col)), drop = FALSE]
  rownames(out) <- NULL
  out
}
