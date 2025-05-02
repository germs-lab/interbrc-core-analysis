brc_report <- function(BRC) {
  quarto::quarto_render(
    input = "brc_report.qmd",
    execute_params = list(BRC = BRC),
    output_format = "html",
    output_file = glue::glue("{BRC}.html")
  )
}
