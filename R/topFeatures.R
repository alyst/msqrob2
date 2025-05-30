#' Toplist of DE proteins, peptides or features
#'
#' @description Summary table of the differentially expressed Features
#'
#' @param models A list with elements of the class `StatModel` that are
#'        estimated using the \code{\link{msqrob}} function
#' @param contrast `numeric` (matrix)vector specifying one contrast of
#'        the linear model coefficients to be tested equal to zero.
#'        The (row)names of the vector should be equal to the names of
#'        parameters of the model.
#' @param adjust.method `character` specifying the method to adjust
#'        the p-values for multiple testing.
#'        Options, in increasing conservatism, include ‘"none"’,
#'        ‘"BH"’, ‘"BY"’ and ‘"holm"’.  See ‘p.adjust’ for the complete
#'        list of options. Default is "BH" the Benjamini-Hochberg method
#'        to controle the False Discovery Rate (FDR).
#' @param sort `boolean(1)` to indicate if the features have to be sorted according
#'        to statistical significance.
#' @param alpha `numeric` specifying the cutoff value for adjusted p-values.
#'        Only features with lower p-values are listed.
#'
#' @examples
#' data(pe)
#'
#' # Aggregate peptide intensities in protein expression values
#' pe <- aggregateFeatures(pe, i = "peptide", fcol = "Proteins", name = "protein")
#'
#' # Fit msqrob model
#' pe <- msqrob(pe, i = "protein", formula = ~condition)
#'
#' # Define contrast
#' getCoef(rowData(pe[["protein"]])$msqrobModels[[1]])
#'
#' # Assess log2 fold change between condition c and condition b:
#' L <- makeContrast("conditionc - conditionb=0", c("conditionb", "conditionc"))
#' topDeProteins <- topFeatures(rowData(pe[["protein"]])$msqrobModels, L)
#' @return A dataframe with log2 fold changes (logFC), standard errors (se),
#'         degrees of freedom of the test (df), t-test statistic (t),
#'         p-values (pval) and adjusted pvalues (adjPval) using the specified
#'         adjust.method in the p.adjust function of the stats package.
#' @importFrom stats pt p.adjust na.exclude
#' @rdname topTable
#' @author Lieven Clement
#' @export

topFeatures <- function(models, contrast, adjust.method = "BH", fix_lmm_ddf = FALSE, sort = TRUE, alpha = 1) {
    if (!is(contrast, "matrix")) {
        contrast <- as.matrix(contrast)
    }
    if (ncol(contrast) > 1) {
        stop("The 'contrast' argument has more than one column, only one contrast is allowed")
    }
    # remove unused coefficients
    contrast <- contrast[rowSums(contrast) != 0, , drop = FALSE]

    logFC <- vapply(models,
        getContrast,
        numeric(1),
        L = contrast
    )
    se <- sqrt(vapply(models,
        varContrast,
        numeric(1),
        L = contrast
    ))
    df <- vapply(models, getDfPosterior, numeric(1))
    ddf_kr <- NULL
    if (fix_lmm_ddf) {
        if (!requireNamespace("parameters", quietly = TRUE)) {
            warning("parameters package is required to calculate the denominator DF for fixed effects in linear mixed models.")
        } else {
            ddf_kr <- bplapply(models, function(model) {
                if ("model" %in% names(model@params) && is(model@params$model, "lmerMod")) {
                    lmm <- model@params$model
                    tryCatch(parameters::dof_kenward(lmm),
                             error = function(e) NA_real_)
                } else {
                    NA_real_
                }
            })
            ddf_kr <- unlist(ddf_kr, recursive = FALSE, use.names = FALSE)
            if (all(is.na(ddf_kr))) {
                ddf_kr <- NULL
            }
        }
    }
    t <- logFC / se
    pval <- pt(-abs(t), df) * 2
    adjPval <- p.adjust(pval, method = adjust.method)
    out <- data.frame(logFC, se, df, t, pval, adjPval)
    if (!is.null(ddf_kr)) {
        out$ddf_kr <- ddf_kr
        out$pval_kr <- pt(-abs(t), ddf_kr) * 2
        out$adjPval_kr <- p.adjust(out$pval_kr, method = adjust.method)
    }
    if (alpha < 1) {
        signif <- adjPval < alpha
        out <- na.exclude(out[signif, ])
    }
    if (sort) {
        out <- out[order(out$pval), ]
    }
    return(out)
}
