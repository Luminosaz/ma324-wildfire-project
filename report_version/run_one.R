#' Run a single AMPL firebreak solve
#'
#' @param model      Character "A" or "B"
#' @param B          Integer budget
#' @param param_name Character "tau" or "alpha"
#' @param param_val  Numeric value for the swept parameter
#'
#' @return A one-row data frame with columns: model, B, param_name,
#'   param_value, damage, n_cleared, cleared_cells, status
run_one = function(model, B, param_name, param_val) {
  mod_file = sprintf("step2_model%s.mod", model)
  run_file = sprintf("step2_model%s.run", model)

  tmp = tempfile(fileext = ".run")
  lines = c(
    "reset;",
    sprintf("model %s;", mod_file),
    "data step2.dat;",
    sprintf("let B := %d;", B),
    sprintf("let %s := %g;", param_name, param_val),
    sprintf("include %s;", run_file)
  )
  writeLines(lines, tmp)

  output = system2(ampl_bin, tmp, stdout = TRUE, stderr = TRUE)
  unlink(tmp)

  # Parse SUMMARY line
  parts = strsplit(grep("^SUMMARY,", output, value = TRUE), ",")[[1]]
  damage = as.numeric(parts[2])
  status = as.integer(parts[3])

  # Parse CLEARED lines
  cleared = sub("^CLEARED,", "", grep("^CLEARED,", output, value = TRUE))

  data.frame(
    model        = model,
    B            = B,
    param_name   = param_name,
    param_value  = param_val,
    damage       = damage,
    n_cleared    = length(cleared),
    cleared_cells = paste(cleared, collapse = ";"),
    status       = status,
    stringsAsFactors = FALSE
  )
}
