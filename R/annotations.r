# ===========================================================================
# R/annotations.R — User-customizable annotation helpers and loader
# ===========================================================================

#' Load Built-in Curated Psychiatric Annotations
#'
#' Retrieves the frozen annotation matrix provided in Supplementary Table S1.21.
#'
#' @param key_type Character, identifier column to use: "entrez" (ID), "ensembl" (ENSEMBL_GENE_ID), or "symbol" (GeneName).
#' @return A data frame containing gene identifiers and 21 biological annotation columns.
#' @export
#' @examples
#' annot <- load_builtin_annotations("entrez")
#' head(annot[, 1:5])
load_builtin_annotations <- function(key_type = c("entrez", "ensembl", "symbol")) {
  key_type <- match.arg(key_type)
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
