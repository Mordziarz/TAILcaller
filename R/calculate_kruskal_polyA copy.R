calculate_polyA_stat_n3 <- function(
  polyA_table     = get_gene_id_out,
  grouping_factor = "group",
  which_level     = "gene_id",
  padj_method     = "fdr"
) {
  if (missing(polyA_table)) {
    stop("'polyA_table' must be defined.")
  }
  required_cols <- c("polyA_length", grouping_factor, which_level)
  if (!all(required_cols %in% colnames(polyA_table))) {
    stop(sprintf(
      "polyA_table must contain columns: %s",
      paste(required_cols, collapse = ", ")
    ))
  }

  message("Starting calculations...")
  start_time <- Sys.time()

  units   <- unique(polyA_table[[which_level]])
  n_units <- length(units)

  out_df <- data.frame(
    test_used = character(n_units),
    p_value   = numeric(n_units),
    stringsAsFactors = FALSE
  )
  out_df[[which_level]] <- units

  for (i in seq_along(units)) {
    unit_id <- units[i]
    subdf   <- polyA_table[polyA_table[[which_level]] == unit_id, , drop = FALSE]
    n_obs   <- nrow(subdf)
    grp_tbl <- table(subdf[[grouping_factor]])
    n_groups <- length(grp_tbl)

    p_val      <- NA_real_
    chosen_test<- NA_character_

    if (n_groups >= 2) {
      if (n_obs < 3) {
        p_val       <- suppressWarnings(
                         stats::kruskal.test(
                           as.formula(paste0("polyA_length ~ ", grouping_factor)),
                           data = subdf
                         )$p.value
                       )
        chosen_test <- "Kruskal–Wallis"
      } else {
        p_norm <- if (n_obs <= 5000) {
          stats::shapiro.test(subdf$polyA_length)$p.value
        } else {
          nortest::lillie.test(subdf$polyA_length)$p.value
        }
        
        if (all(grp_tbl >= 2)) {
          lev_p <- tryCatch({
            car::leveneTest(
              as.formula(paste0("polyA_length ~ ", grouping_factor)),
              data = subdf
            )[["Pr(>F)"]][1]
          }, error = function(e) NA_real_)
        } else {
          lev_p <- NA_real_
        }

        if (n_groups > 2 && !is.na(p_norm) && p_norm > 0.05) {
          # Dane są normalne
          if (!is.na(lev_p) && lev_p > 0.05) {
            aov_res <- stats::aov(
              as.formula(paste0("polyA_length ~ ", grouping_factor)),
              data = subdf
            )
            p_val       <- summary(aov_res)[[1]][["Pr(>F)"]][1]
            chosen_test <- "ANOVA"
          } else {
            welch_res <- stats::oneway.test(
              as.formula(paste0("polyA_length ~ ", grouping_factor)),
              data = subdf,
              var.equal = FALSE
            )
            p_val       <- welch_res$p.value
            chosen_test <- "Welch ANOVA"
          }
        } else {
          p_val       <- suppressWarnings(
                           stats::kruskal.test(
                             as.formula(paste0("polyA_length ~ ", grouping_factor)),
                             data = subdf
                           )$p.value
                         )
          chosen_test <- "Kruskal–Wallis"
        }
      }
    }

    out_df$test_used[i] <- chosen_test
    out_df$p_value[i]   <- p_val
  }

  out_df$padj <- stats::p.adjust(out_df$p_value, method = padj_method)

  end_time <- Sys.time()
  message(
    "Done in ",
    round(difftime(end_time, start_time, units = "mins"), 2),
    " mins."
  )

  out_df <- out_df[c(which_level, "test_used", "p_value", "padj")]
  return(out_df)
}
