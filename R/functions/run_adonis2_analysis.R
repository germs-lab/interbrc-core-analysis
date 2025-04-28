run_adonis2_analysis <- function(
  ps_object,
  formula,
  strata = NULL,
  nperm = 1000,
  method = "bray"
) {
  # Create relative abundance
  ps_rel <- transform_sample_counts(ps_object, function(x) x / sum(x))

  # Calculate Bray-Curtis distances
  BC.dist <- vegan::vegdist(t(otu_table(ps_rel)), method = method)

  # Sample data for models
  sample_data_df <- data.frame(sample_data(ps_rel))

  # Handle the formula properly
  if (is.character(formula)) {
    formula <- as.formula(paste("BC.dist ~", formula))
  } else {
    # If it's already a formula, update it to include BC.dist as response
    formula_text <- paste(
      "BC.dist ~",
      paste(as.character(formula)[-1], collapse = " ")
    )
    formula <- as.formula(formula_text)
  }

  # Run adonis2 with the provided formula
  if (is.null(strata)) {
    # Without strata
    result <- vegan::adonis2(
      formula,
      data = sample_data_df,
      nperm = nperm,
      method = method,
      by = "terms"
    )

    return(list(full_model = result))
  } else {
    # With strata
    strata_var <- sample_data_df[[strata]]

    # Run model without strata
    result <- vegan::adonis2(
      formula,
      data = sample_data_df,
      nperm = nperm,
      method = method,
      by = "terms"
    )

    # Run model with strata
    result2 <- vegan::adonis2(
      formula,
      data = sample_data_df,
      nperm = nperm,
      method = method,
      by = "terms",
      strata = strata_var
    )

    return(list(
      full_model = result,
      treatment_model_strata = result2
    ))
  }
}
