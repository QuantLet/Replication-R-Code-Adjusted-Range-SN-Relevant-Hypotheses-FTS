rm(list = ls())

# Build the CPP change-point estimator histograms (fig:CPP_estimate_change_point)
# Requires simulation tables from output/simulations/CPP/fIID (produced by run_simulations_CPP_fIID.R)

library(fda)
source(file.path("R", "generalMethods.R"))

aVector <- 0:25/50
ntables <- length(aVector)

base_sim <- file.path("output", "simulations", "CPP", "fIID", "n200")
out_dir <- file.path("output", "figures", "CPP")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# Read baseline tables (column 3 = change-point estimator)
n200list <- readTables(base_sim, ntables)

# Indices for a = 0.1, 0.2, 0.3 given aVector = 0, 0.02, ..., 0.5
idx_a01 <- which(abs(aVector - 0.1) < 1e-12)
idx_a02 <- which(abs(aVector - 0.2) < 1e-12)
idx_a03 <- which(abs(aVector - 0.3) < 1e-12)

plot_hist <- function(cp_vec, filename) {
  png(file.path(out_dir, filename),
      width = 5.65, height = 5.25, units = "cm", res = 1200, pointsize = 4)
  histo <- hist(cp_vec, main = "", xlab = "", breaks = 1:100 * 2, font.lab = 2, lwd = 0.5)
  par(lwd = 0.5)
  plot(histo, bty = "n", axes = FALSE, main = "", xlab = "")
  box(lwd = 0.5)
  axis(side = 1, lwd = 0.5)
  axis(side = 2, lwd = 0.5)
  dev.off()
}

plot_hist(n200list[[idx_a01]][, 3], "a01CPest.png")
plot_hist(n200list[[idx_a02]][, 3], "a02CPest.png")
plot_hist(n200list[[idx_a03]][, 3], "a03CPest.png")

message("CPP change-point estimator histograms written to: ", normalizePath(out_dir))
