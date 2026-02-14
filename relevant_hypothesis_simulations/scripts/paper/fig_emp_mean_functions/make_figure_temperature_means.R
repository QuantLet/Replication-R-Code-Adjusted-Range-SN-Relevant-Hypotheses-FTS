rm(list = ls())

# Reproduce Figure fig:emp-mean-functions (mean annual temperature curves)
# Assumes working directory is the repository root.

library(fda)
library(rnoaa)
library(dplyr)
library(lubridate)

# -----------------------------
# Stations and period (as in manuscript)
# -----------------------------
station_portland <- "USW00014764"   # Portland Jetport
station_seattle  <- "USW00024233"   # Seattle-Tacoma International Airport

date_min <- "1940-11-01"
date_max <- "2024-12-31"

nbasis <- 49
day_grid <- 1:365

# -----------------------------
# Helpers
# -----------------------------
pull_tmin <- function(station_id, date_min, date_max) {
  # rnoaa returns temps in tenths of °C
  df <- meteo_pull_monitors(monitors = station_id, date_min = date_min, date_max = date_max, var = "TMIN") %>%
    transmute(
      date = as.Date(date),
      tmin = tmin / 10
    )
  df
}

build_year_matrix <- function(df) {
  df <- df %>%
    mutate(
      year = year(date),
      month = month(date),
      day = day(date),
      leap = leap_year(date),
      doy = yday(date)
    ) %>%
    # drop Feb 29 so all years have 365 days
    filter(!(month == 2 & day == 29)) %>%
    mutate(
      # in leap years, shift days after Feb 29 by -1
      doy365 = ifelse(leap & month > 2, doy - 1, doy)
    )

  years <- sort(unique(df$year))
  Y <- matrix(NA_real_, nrow = 365, ncol = length(years))
  colnames(Y) <- years

  for (j in seq_along(years)) {
    yr <- years[j]
    sub <- df %>% filter(year == yr)

    v <- rep(NA_real_, 365)
    v[sub$doy365] <- sub$tmin

    # simple interpolation for missing days
    if (anyNA(v)) {
      idx <- which(!is.na(v))
      if (length(idx) >= 2) {
        v <- approx(x = idx, y = v[idx], xout = 1:365, rule = 2)$y
      }
    }
    Y[, j] <- v
  }

  list(Y = Y, years = years)
}

mean_fd_from_matrix <- function(Y, nbasis) {
  basis <- create.fourier.basis(rangeval = c(1, 365), nbasis = nbasis)
  fdobj <- Data2fd(argvals = day_grid, y = Y, basisobj = basis)
  mean_coef <- rowMeans(fdobj$coefs)
  fd(mean_coef, basis)
}

# -----------------------------
# Download and build FD objects
# -----------------------------
message("Downloading Portland data ...")
port_df <- pull_tmin(station_portland, date_min, date_max)

message("Downloading Seattle data ...")
sea_df <- pull_tmin(station_seattle, date_min, date_max)

port_mat <- build_year_matrix(port_df)
sea_mat  <- build_year_matrix(sea_df)

port_mean_fd <- mean_fd_from_matrix(port_mat$Y, nbasis = nbasis)
sea_mean_fd  <- mean_fd_from_matrix(sea_mat$Y, nbasis = nbasis)

port_mean <- as.numeric(eval.fd(day_grid, port_mean_fd))
sea_mean  <- as.numeric(eval.fd(day_grid, sea_mean_fd))

# -----------------------------
# Plot
# -----------------------------
out_dir <- file.path("output", "figures")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

png(file.path(out_dir, "emp_mean_functions.png"),
    width = 9, height = 6, units = "cm", res = 600, pointsize = 8)

plot(day_grid, port_mean, type = "l",
     xlab = "Day of year", ylab = "Mean daily minimum temperature (°C)",
     lwd = 1, bty = "n", axes = FALSE)

box(lwd = 0.5)
axis(side = 1, lwd = 0.5)
axis(side = 2, lwd = 0.5)

lines(day_grid, sea_mean, lty = 2, lwd = 1)

legend("topleft",
       legend = c("Portland Jetport", "Seattle–Tacoma Int'l Airport"),
       lty = c(1, 2), lwd = 1, bty = "n")

dev.off()

message("Figure written to: ", normalizePath(file.path(out_dir, "emp_mean_functions.png")))
