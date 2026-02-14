rm(list = ls())

# Reproduce Figure fig:gbp-comparison (GBP/USD mean curves before/after Brexit day)
#
# IMPORTANT:
# This script does NOT ship data. You must provide your own CSV file.
#
# Expected columns (minimum):
#   - DateTime : timestamp (e.g. "2016-01-04 00:01:00")
#   - Close    : GBP/USD close price (USD per GBP)
#
# Output:
#   - output/figures/gbp-usd-before.png
#   - output/figures/gbp-usd-after.png

library(fda)
library(dplyr)
library(lubridate)

# -----------------------------
# User input
# -----------------------------
csv_path <- "PATH/TO/YOUR/GBPUSD_2016.csv"   # <-- EDIT THIS
brexit_date <- as.Date("2016-06-23")

# -----------------------------
# Load data
# -----------------------------
if (!file.exists(csv_path)) {
  stop("CSV file not found. Please edit csv_path in this script.\n  csv_path = ", csv_path)
}

dat <- read.csv(csv_path, stringsAsFactors = FALSE)

# Basic column checks
if (!("Close" %in% names(dat))) stop("CSV must contain a 'Close' column.")
if (!("DateTime" %in% names(dat))) stop("CSV must contain a 'DateTime' column (timestamp).")

dat <- dat %>%
  mutate(
    DateTime = ymd_hms(DateTime, tz = "UTC"),
    Date = as.Date(DateTime),
    # time of day scaled to [0, 1]
    DayTime = (hour(DateTime) * 3600 + minute(DateTime) * 60 + second(DateTime)) / (24 * 3600),
    # convert to USD/GBP if needed (original code used 1/Close * 100 for scaling)
    USD_GBP = (1 / Close) * 100
  ) %>%
  arrange(Date, DayTime)

# -----------------------------
# Build daily functional data (B-splines)
# -----------------------------
daily_list <- split(dat, dat$Date)

daily_counts <- sapply(daily_list, nrow)
min_obs <- min(daily_counts)

# Heuristic basis size (as in the original project script)
nbasis <- min(21, floor(min_obs * 0.8))
message("Minimum observations per day: ", min_obs)
message("Using nbasis = ", nbasis)

basis <- create.bspline.basis(rangeval = c(0, 1), nbasis = nbasis)

n_days <- length(daily_list)
coefs <- matrix(NA_real_, nrow = nbasis, ncol = n_days)
day_names <- names(daily_list)

for (i in seq_along(daily_list)) {
  day_df <- daily_list[[i]]

  # Smooth within the day (irregular timestamps allowed)
  fit <- smooth.basis(argvals = day_df$DayTime, y = day_df$USD_GBP, fdParobj = basis)
  coefs[, i] <- fit$fd$coefs
}

fdata <- fd(coefs, basis)
fdata$fdnames$reps <- day_names

# Identify Brexit index
date_vec <- as.Date(day_names)
brexit_idx <- which(date_vec == brexit_date)
if (length(brexit_idx) != 1) {
  stop("Could not find exactly one Brexit date (", brexit_date, ") in the data. Found: ", length(brexit_idx))
}

mean_before <- mean.fd(fdata[1:(brexit_idx - 1)])
mean_after  <- mean.fd(fdata[brexit_idx:length(date_vec)])

# Evaluate mean curves
points <- 1:101 / 101
mean_before_vals <- as.numeric(eval.fd(points, mean_before))
mean_after_vals  <- as.numeric(eval.fd(points, mean_after))

# -----------------------------
# Plot
# -----------------------------
out_dir <- file.path("output", "figures")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# Before
png(file.path(out_dir, "gbp-usd-before.png"),
    width = 5.65, height = 5.25, units = "cm", res = 1200, pointsize = 4)

plot(points, mean_before_vals, type = "l",
     xlab = expression(t), ylab = "USD / GBP (scaled)",
     font.lab = 2, lwd = 0.5, cex.lab = 1, bty = "n", axes = FALSE)

box(lwd = 0.5)
axis(side = 1, lwd = 0.5)
axis(side = 2, lwd = 0.5)
dev.off()

# After
png(file.path(out_dir, "gbp-usd-after.png"),
    width = 5.65, height = 5.25, units = "cm", res = 1200, pointsize = 4)

plot(points, mean_after_vals, type = "l",
     xlab = expression(t), ylab = "USD / GBP (scaled)",
     font.lab = 2, lwd = 0.5, cex.lab = 1, bty = "n", axes = FALSE)

box(lwd = 0.5)
axis(side = 1, lwd = 0.5)
axis(side = 2, lwd = 0.5)
dev.off()

message("GBP/USD figures written to: ", normalizePath(out_dir))
