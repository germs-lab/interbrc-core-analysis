core_select_pipeline <- function(physeq_obj, .title, .subtitle) {
  # Core selection
  core_summary_lists <- extract_core(
    physeq_obj,
    Var = "site",
    method = "increase",
    increase_value = 2
  )

  # Plot Bray-Curtis Dissimilarity Curve:
  bc_curve <- brc_bc_curve(core_summary_list = core_summary_lists)

  occabun_plot <- brc_bc_occ_curve(
    core_summary_list = core_summary_lists
  )

  # Print
  bc_curve +
    ggtitle(
      .title,
      .subtitle
    )

  occabun_plot +
    ggtitle(
      .title,
      .subtitle
    )

  return(list(
    core_summary = core_summary_lists,
    bc_curve_threshold = bc_curve,
    occupancy_plot = occabun_plot
  ))
}
