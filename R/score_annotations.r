# R/score_annotations.r — Automated scoring for DAVID biological annotations
# ===========================================================================

#' Default Keyword Dictionaries for Biological Annotation Scoring
#'
#' Provides the deterministic keyword lists documented in Supplementary Methods Section 1.
#'
#' @return A named list of character vectors for each annotation category.
#' @export
#' @examples
#' kw <- default_annotation_keywords()
#' names(kw)
default_annotation_keywords <- function() {
  list(
    expre_brain = c(
      "brain", "cortex", "cerebellar", "cerebellum", "putamen", "caudate",
      "nucleus accumbens", "hippocampus", "amygdala", "spinal cord",
      "hypothalamus", "substantia nigra", "synapse", "synaptic"
    ),
    dise_brain = c(
      "schizo", "depress", "bipolar", "autis", "attention deficit",
      "alcohol", "circadian", "drug abuse"
    ),
    pathway_brain = c(
      "neuro", "alcohol", "synapse", "synaptic", "serotonergic", "gabaergic"
    ),
    bio_brain = c(
      "brain", "dopamine", "serotonin", "neuro", "tryptophan", "epinephrine",
      "histamine", "norepinephrine", "circadian", "cortex", "cerebellar",
      "cerebellum", "putamen", "caudate", "nucleus accumbens", "hippocampus",
      "amygdala", "spinal cord", "hypothalamus", "substantia nigra", "schizo",
      "depress", "bipolar", "autis", "attention deficit", "alcohol",
      "synapse", "synaptic"
    ),
    bio_inter = c(
      "brain", "dopamine", "serotonin", "neuro", "tryptophan", "epinephrine",
      "histamine", "norepinephrine", "circadian", "synapse", "synaptic"
    )
  )
}

#' Helper: Deterministic binary scoring rule (0 or 1)
#' @noRd
.score_binary_text <- function(text_vec, keywords) {
  text_vec <- as.character(text_vec)
  text_vec[is.na(text_vec)] <- ""
  text_lower <- tolower(trimws(text_vec))
  
  # Build regex pattern for substring match
  escaped_kw <- paste0(keywords, collapse = "|")
  has_match <- grepl(escaped_kw, text_lower, perl = TRUE)
  as.integer(has_match & nzchar(text_lower))
}

#' Helper: Deterministic tiered scoring rule (0, 1, or 2)
#' @noRd
.score_tiered_text <- function(text_vec, core_keywords) {
  text_vec <- as.character(text_vec)
  text_vec[is.na(text_vec)] <- ""
  text_lower <- tolower(trimws(text_vec))
  
  res <- integer(length(text_vec)) # defaults to 0
  is_annotated <- nzchar(text_lower) & !text_lower %in% c("na", "null", "none", "-", ".")
  res[is_annotated] <- 1L
  
  escaped_kw <- paste0(core_keywords, collapse = "|")
  has_core <- grepl(escaped_kw, text_lower, perl = TRUE)
  res[is_annotated & has_core] <- 2L
  res
}

#' Helper: TFBS count / text scoring rule (0 or 1)
#' @noRd
.score_tfbs_val <- function(tfbs_vec) {
  res <- integer(length(tfbs_vec))
  for (i in seq_along(tfbs_vec)) {
    val <- tfbs_vec[i]
    if (is.na(val) || is.null(val)) next
    
    num_val <- suppressWarnings(as.numeric(val))
    if (!is.na(num_val)) {
      if (num_val > 1) res[i] <- 1L
    } else {
      # Parse text list separated by commas or semicolons
      s <- trimws(as.character(val))
      if (nzchar(s) && !tolower(s) %in% c("na", "null", "none", "-", ".")) {
        parts <- unlist(strsplit(gsub(";", ",", s), ","))
        parts <- trimws(parts)
        parts <- parts[nzchar(parts)]
        if (length(parts) > 1) res[i] <- 1L
      }
    }
  }
  res
}

#' Compute Deterministic Biological Annotation Scores from Raw Database Outputs
#'
#' Transforms raw text/term query outputs (e.g. from DAVID 6.8 or similar functional databases)
#' into the 6 deterministic tiered and binary numerical features (`expre.brain`, `dise.brain`,
#' `pathway.brain`, `tfbs`, `bio.brain`, `bio.inter`) as specified in Supplementary Methods Section 1.
#'
#' @param raw_df A data frame containing raw annotation columns and an identifier column.
#' @param id_col Character, name of the gene identifier column (e.g. "ID", "GENE", "gene").
#' @param col_mapping Optional named list mapping standard feature names to column names in `raw_df`.
#'   Standard names: \code{tissue} (for \code{expre.brain}), \code{disease} (for \code{dise.brain}),
#'   \code{pathway} (for \code{pathway.brain}), \code{tfbs} (for \code{tfbs}), \code{process}
#'   (for \code{bio.brain}), \code{interaction} (for \code{bio.inter}).
#' @param keywords Optional custom keyword dictionary list overriding \code{\link{default_annotation_keywords}}.
#' @return A data frame containing the \code{id_col} and the scored numeric columns ready for FDRreg.
#' @export
#' @examples
#' raw_data <- data.frame(
#'   ID = c("1001", "1002", "1003"),
#'   tissue = c("Fetal brain, cortex", "Liver", NA),
#'   disease = c("Schizophrenia, Bipolar", "Diabetes", NA),
#'   pathway = c("Dopaminergic synapse", "Glycolysis", NA),
#'   tfbs_count = c(4, 0, NA),
#'   process = c("Synaptic transmission", "Cell cycle", NA),
#'   interaction = c("Dopamine receptor binding", "Insulin binding", NA)
#' )
#' scores <- score_david_annotations(raw_data, id_col = "ID")
#' head(scores)
score_david_annotations <- function(raw_df, id_col = "ID",
                                    col_mapping = list(
                                      tissue = "tissue",
                                      disease = "disease",
                                      pathway = "pathway",
                                      tfbs = "tfbs_count",
                                      process = "process",
                                      interaction = "interaction"
                                    ),
                                    keywords = NULL) {
  if (!is.data.frame(raw_df)) {
    stop("'raw_df' must be a data frame.", call. = FALSE)
  }
  if (!id_col %in% colnames(raw_df)) {
    stop(sprintf("Identifier column '%s' not found in input data frame.", id_col), call. = FALSE)
  }
  
  kw <- if (is.null(keywords)) default_annotation_keywords() else keywords
  
  res <- data.frame(row.names = seq_len(nrow(raw_df)))
  res[[id_col]] <- as.character(raw_df[[id_col]])
  
  # 1. expre.brain (binary)
  col_tissue <- col_mapping$tissue
  if (!is.null(col_tissue) && col_tissue %in% colnames(raw_df)) {
    res$expre.brain <- .score_binary_text(raw_df[[col_tissue]], kw$expre_brain)
  }
  
  # 2. dise.brain (tiered 0, 1, 2)
  col_dise <- col_mapping$disease
  if (!is.null(col_dise) && col_dise %in% colnames(raw_df)) {
    res$dise.brain <- .score_tiered_text(raw_df[[col_dise]], kw$dise_brain)
  }
  
  # 3. pathway.brain (tiered 0, 1, 2)
  col_path <- col_mapping$pathway
  if (!is.null(col_path) && col_path %in% colnames(raw_df)) {
    res$pathway.brain <- .score_tiered_text(raw_df[[col_path]], kw$pathway_brain)
  }
  
  # 4. tfbs (binary)
  col_tfbs <- col_mapping$tfbs
  if (!is.null(col_tfbs) && col_tfbs %in% colnames(raw_df)) {
    res$tfbs <- .score_tfbs_val(raw_df[[col_tfbs]])
  }
  
  # 5. bio.brain (tiered 0, 1, 2)
  col_proc <- col_mapping$process
  if (!is.null(col_proc) && col_proc %in% colnames(raw_df)) {
    res$bio.brain <- .score_tiered_text(raw_df[[col_proc]], kw$bio_brain)
  }
  
  # 6. bio.inter (tiered 0, 1, 2)
  col_inter <- col_mapping$interaction
  if (!is.null(col_inter) && col_inter %in% colnames(raw_df)) {
    res$bio.inter <- .score_tiered_text(raw_df[[col_inter]], kw$bio_inter)
  }
  
  res
}
