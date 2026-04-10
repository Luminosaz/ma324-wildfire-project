# Project R profile — auto-loaded when R starts in this directory.

# Make AMPL binary findable when calling via system() / system2().
Sys.setenv(
  AMPL_BIN    = "/Users/jerrymac/Library/Python/3.9/lib/python/site-packages/ampl_module_base/bin/ampl",
  GUROBI_BIN  = "/Users/jerrymac/Library/Python/3.9/lib/python/site-packages/ampl_module_gurobi/bin/gurobi",
  HIGHS_BIN   = "/Users/jerrymac/Library/Python/3.9/lib/python/site-packages/ampl_module_highs/bin/highs",
  CBC_BIN     = "/Users/jerrymac/Library/Python/3.9/lib/python/site-packages/ampl_module_cbc/bin/cbc"
)

# Add AMPL to PATH so `ampl` works in shells spawned by R.
old_path <- Sys.getenv("PATH")
ampl_dir <- "/Users/jerrymac/Library/Python/3.9/lib/python/site-packages/ampl_module_base/bin"
if (!grepl(ampl_dir, old_path, fixed = TRUE)) {
  Sys.setenv(PATH = paste(ampl_dir, old_path, sep = ":"))
}

cat("MA324 project R profile loaded.\n")
cat("AMPL_BIN   =", Sys.getenv("AMPL_BIN"), "\n")
cat("GUROBI_BIN =", Sys.getenv("GUROBI_BIN"), "\n")
cat("HIGHS_BIN  =", Sys.getenv("HIGHS_BIN"), "\n")
