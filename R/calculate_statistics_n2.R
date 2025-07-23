#' Calculate Statistical Differences Between Two Groups for Each Molecule
#'
#' For each unique molecule (defined by `which_level`), this function compares
#' poly(A) tail lengths between a `control_group` and a `treated_group` using
#' appropriate statistical tests:
#' \itemize{
#'   \item \strong{Student's t-test}: If data in both groups exhibit a normal distribution and have homogeneous variances.
#'   \item \strong{Welch's t-test}: If data in both groups exhibit a normal distribution but have unequal variances.
#'   \item \strong{Wilcoxon rank-sum test (Mann-Whitney U test)}: Otherwise (i.e., non-normality, or fewer than 2 observations in any group).
#' }
#' It also calculates mean tail lengths for each group, fold change,
#' difference in length, Cohen's d, and adjusts p-values for multiple comparisons.
#'
#' @param polyA_table A \code{data.frame} with per-read poly(A) tail lengths. Must contain:
#'   \describe{
#'     \item{\code{polyA_length} (numeric)}{Tail length of each read.}
#'     \item{\code{<grouping_factor>} (character or factor)}{Grouping variable for comparison, e.g., experimental condition.}
#'     \item{\code{<which_level>} (character or factor)}{Molecular identifier, e.g., gene_id, transcript_id.}
#'   }
#' @param grouping_factor Character; name of the column in \code{polyA_table} to use for defining groups. Defaults to \code{"group"}.
#' @param which_level Character; name of the column in \code{polyA_table} that uniquely identifies each molecule for which statistics should be calculated. Defaults to \code{"gene_id"}.
#' @param control_group Character; name of the control group level in `grouping_factor`.
#' @param treated_group Character; name of the treated group level in `grouping_factor`.
#' @param padj_method Character; method for p-value adjustment. See `p.adjust.methods` for options. Defaults to \code{"fdr"}.
#'
#' @return A \code{data.frame} with one row per unique molecule, containing:
#'   \describe{
#'     \item{\code{<which_level>} (character/factor)}{Molecule identifier.}
#'     \item{\code{p_value} (numeric)}{P-value from the chosen statistical test.}
#'     \item{\code{mean_group_ctr} (numeric)}{Mean poly(A) length in the control group.}
#'     \item{\code{mean_group_trt} (numeric)}{Mean poly(A) length in the treated group.}
#'     \item{\code{diff_length} (numeric)}{Difference in mean poly(A) length (treated - control).}
#'     \item{\code{fold_change} (numeric)}{Fold change of mean poly(A) length (treated / control).}
#'     \item{\code{cohen_d} (numeric)}{Cohen's d effect size.}
#'     \item{\code{cohen_effect} (character)}{Categorical interpretation of Cohen's d ("small", "medium", "large").}
#'     \item{\code{padj} (numeric)}{Adjusted p-value.}
#'     \item{\code{Log2FC} (numeric)}{Log2 Fold Change.}
#'     \item{\code{test_performed} (character)}{Name of the statistical test performed for the molecule.}
#'   }
#'
#' @details
#' This function iterates through each unique molecule identified by `which_level`.
#' For each molecule, it extracts poly(A) tail lengths for the specified control and
#' treated groups. It then applies the following logic for statistical testing:
#' \enumerate{
#'   \item \strong{Minimum Observations Check}: If either the control or treated group has fewer than 2 observations for a given molecule, `NA`s are returned for statistical results, as a standard deviation cannot be calculated for a single observation, making t-tests or Cohen's d impossible.
#'   \item \strong{Normality Test}: Shapiro-Wilk test (for N <= 5000) or Lilliefors test (for N > 5000) is performed on both control and treated data subsets for the current molecule.
#'   \item \strong{Homogeneity of Variance Test}: Levene's test is performed if both groups have at least 2 observations.
#'   \item \strong{Test Selection}:
#'     \itemize{
#'       \item If *both* groups are normally distributed *and* have homogeneous variances, a \strong{Student's t-test} (`var.equal = TRUE`) is conducted.
#'       \item If *both* groups are normally distributed *but* have unequal variances, a \strong{Welch's t-test} (`var.equal = FALSE`) is conducted.
#'       \item Otherwise (if either group is non-normal, or fewer than 3 observations for normality test, or Levene's test cannot be performed, or variances are unequal and normality doesn't hold for both), a \strong{Wilcoxon rank-sum test} is conducted.
#'     \item \strong{Effect Size}: Cohen's d is calculated based on the pooled standard deviation (for cases where it's applicable, primarily for parametric tests, but included for all as a general effect size measure).
#'     \item \strong{Multiple Comparison Adjustment}: P-values are adjusted using the specified `padj_method`.
#'   \end{enumerate}
#'
#' @author Mateusz Mazdziarz
#'
#' @importFrom stats wilcox.test t.test shapiro.test p.adjust sd mean
#' @importFrom nortest lillie.test
#' @importFrom car leveneTest
#' @importFrom dplyr %>% distinct
#' @importFrom rlang sym
#' @export

calculate_statistics_n2 <- function(polyA_table = get_gene_id_out, grouping_factor = "group", which_level = "gene_id", control_group = NULL, treated_group = NULL, padj_method = "fdr") {
  if (missing(polyA_table)) stop("'polyA_table' must be defined.")
  if (is.null(control_group) | is.null(treated_group)) {
    stop("Both 'control_group' and 'treated_group' must be provided.")
  }

  message("Starting processing...")
  start_time <- Sys.time()

  if (!all(c("polyA_length", grouping_factor, which_level) %in% colnames(polyA_table))) {
    stop(paste0("The 'polyA_table' must contain 'polyA_length', '", grouping_factor, "', and '", which_level, "' columns."))
  }

  polyA_table_filtered <- polyA_table[polyA_table[[grouping_factor]] %in% c(control_group, treated_group), ]

  how_molecules <- polyA_table_filtered %>%
    dplyr::distinct(!!rlang::sym(which_level)) %>%
    as.data.frame()

  results <- lapply(1:nrow(how_molecules), function(i) {
    molecule <- how_molecules[[which_level]][i]
    subset_data <- polyA_table_filtered[polyA_table_filtered[[which_level]] == molecule, ]

    ctr_data <- subset_data$polyA_length[subset_data[[grouping_factor]] == control_group]
    trt_data <- subset_data$polyA_length[subset_data[[grouping_factor]] == treated_group]

    n_ctr <- length(ctr_data)
    n_trt <- length(trt_data)

    test_name <- NA

    if (n_ctr < 2 || n_trt < 2) {
      return(list(
        p_value = NA,
        mean_ctr = NA,
        mean_trt = NA,
        fold_change = NA,
        cohen_d = NA,
        cohen_effect = NA,
        diff_length = NA,
        test_performed = "Insufficient data (<2 observations per group)"
      ))
    }

    shapiro_or_lillie <- function(x, group_name = "") {
      if (length(x) < 3) {
        warning(paste0("Not enough observations (n=", length(x), ") in group '", group_name, "' for normality test. Normality cannot be assessed."))
        return(NA)
      }
      if (length(unique(x)) == 1) {
        warning(paste0("All 'polyA_length' values are identical in group '", group_name, "'. Normality test cannot be performed. Normality cannot be assessed."))
        return(NA)
      }

      tryCatch({
        if (length(x) <= 5000) {
          stats::shapiro.test(x)$p.value
        } else {
          if (!requireNamespace("nortest", quietly = TRUE)) {
            stop("Package 'nortest' is required for Lilliefor's test when n > 5000. Please install it.")
          }
          nortest::lillie.test(x)$p.value
        }
      }, error = function(e) {
        warning(paste0("Error in normality test for group '", group_name, "' (n=", length(x), "): ", e$message, ". Normality cannot be assessed."))
        return(NA)
      })
    }


    normality_ctr_p <- shapiro_or_lillie(ctr_data)
    normality_trt_p <- shapiro_or_lillie(trt_data)

    all_normal <- !is.na(normality_ctr_p) && normality_ctr_p > 0.05 &&
                  !is.na(normality_trt_p) && normality_trt_p > 0.05


    levene_p <- NA
    if (n_ctr >= 2 && n_trt >= 2) {
      temp_data <- data.frame(
        polyA_length = c(ctr_data, trt_data),
        group = factor(c(rep(control_group, n_ctr), rep(treated_group, n_trt)))
      )
      levene_res <- tryCatch(
        {

          car::leveneTest(as.formula(paste0("polyA_length ~ ", grouping_factor)), data = temp_data, center = median)
        },
        error = function(e) {
          warning(paste0("Levene's test failed for molecule '", molecule, "': ", e$message))
          NULL
        }
      )
      if (!is.null(levene_res)) {
        levene_p <- levene_res[["Pr(>F)"]][1]
      }
    }
    homogeneous_variances <- !is.na(levene_p) && levene_p > 0.05

p_val <- NA
    if (all_normal && homogeneous_variances) {
      p_val <- suppressWarnings(
        stats::t.test(ctr_data, trt_data, var.equal = TRUE)$p.value
      )
      test_name <- "Student's t-test"
    } else if (all_normal && !homogeneous_variances) {
      p_val <- suppressWarnings(
        stats::t.test(ctr_data, trt_data, var.equal = FALSE)$p.value
      )
      test_name <- "Welch's t-test"
    } else {
      p_val <- suppressWarnings(
        stats::wilcox.test(ctr_data, trt_data)$p.value
      )
      test_name <- "Wilcoxon rank-sum test"
    }

    mean_ctr <- mean(ctr_data, na.rm = TRUE)
    mean_trt <- mean(trt_data, na.rm = TRUE)

    fc <- ifelse(mean_ctr == 0, NA, mean_trt / mean_ctr)

    sd_ctr <- sd(ctr_data, na.rm = TRUE)
    sd_trt <- sd(trt_data, na.rm = TRUE)

    pooled_sd <- NA
    if (n_ctr >= 2 && n_trt >= 2) {
      pooled_sd <- sqrt(((n_ctr - 1) * sd_ctr^2 + (n_trt - 1) * sd_trt^2) / (n_ctr + n_trt - 2))
    }

    cd <- NA
    if (!is.na(pooled_sd) && pooled_sd > 0) {
      cd <- abs((mean_trt - mean_ctr) / pooled_sd)
    }

    list(
      p_value = p_val,
      mean_ctr = mean_ctr,
      mean_trt = mean_trt,
      fold_change = fc,
      cohen_d = cd,
      cohen_effect = ifelse(is.na(cd), NA,
                            ifelse(cd < 0.2, "small",
                                   ifelse(cd < 0.5, "medium", "large"))),
      diff_length = mean_trt - mean_ctr,
      test_performed = test_name
    )
  })


  how_molecules$p_value <- sapply(results, `[[`, "p_value")
  how_molecules$mean_group_ctr <- sapply(results, `[[`, "mean_ctr")
  how_molecules$mean_group_trt <- sapply(results, `[[`, "mean_trt")
  how_molecules$diff_length <- sapply(results, `[[`, "diff_length")
  how_molecules$fold_change <- sapply(results, `[[`, "fold_change")
  how_molecules$cohen_d <- sapply(results, `[[`, "cohen_d")
  how_molecules$cohen_effect <- sapply(results, `[[`, "cohen_effect")
  how_molecules$test_performed <- sapply(results, `[[`, "test_performed")

  valid_p_values_idx <- !is.na(how_molecules$p_value)
  how_molecules$padj <- NA
  how_molecules$padj[valid_p_values_idx] <- p.adjust(how_molecules$p_value[valid_p_values_idx], method = padj_method)

  how_molecules$Log2FC <- log2(how_molecules$fold_change)

  message("Processing completed. Time: ",
    round(difftime(Sys.time(), start_time, units = "mins"), 1),
    " minutes"
  )

  return(how_molecules)
}