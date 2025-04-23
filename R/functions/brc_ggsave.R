brc_ggsave <- function(objects, save_path) {
  ggnames <- c(names(objects))

  plot_paths <- str_glue("{save_path}{tolower(ggnames)}.png")

  purrr::walk2(
    plot_paths,
    objects,
    \(path, plot)
      ggsave(
        filename = path,
        plot = plot,
        dpi = 300,
        width = 200,
        height = 200,
        units = "mm"
      )
  )
}
