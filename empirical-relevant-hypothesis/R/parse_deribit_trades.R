# Parse Deribit-style BTC options trade CSV (s,t,trade)
#
# The CSV you provided has:
#   - t: millisecond timestamp (also duplicated inside the JSON)
#   - trade: JSON string with fields like instrument_name, iv, amount, index_price
#
# This loader:
#   - reads the CSV with data.table::fread
#   - extracts only the fields we need using fast regex (stringi)
#   - parses instrument_name like "BTC-30JAN26-80000-P" into expiry/strike/type

parse_instrument_name <- function(instr) {
  # Expected: BTC-30JAN26-80000-P
  #           BTC-16JAN26-85000-C
  # Returns data.table with columns: expiry_date, strike, cp_flag

  # Capture groups: day, mon, year, strike, flag
  m <- stringi::stri_match_first_regex(
    instr,
    "^BTC-([0-9]{1,2})([A-Z]{3})([0-9]{2})-([0-9]+(?:\\.[0-9]+)?)-([CP])$"
  )
  day   <- as.integer(m[, 2])
  mon_s <- m[, 3]
  yy    <- as.integer(m[, 4])
  strike <- as.numeric(m[, 5])
  cp    <- m[, 6]

  mon_map <- c(JAN=1, FEB=2, MAR=3, APR=4, MAY=5, JUN=6,
               JUL=7, AUG=8, SEP=9, OCT=10, NOV=11, DEC=12)
  mon <- unname(mon_map[mon_s])
  # assume 20xx
  year <- 2000L + yy
  expiry_date <- as.Date(sprintf("%04d-%02d-%02d", year, mon, day), tz = "UTC")

  data.table::data.table(
    expiry_date = expiry_date,
    strike = strike,
    cp_flag = cp
  )
}


load_deribit_trade_csv <- function(path,
                                   max_rows = Inf,
                                   keep_cols = c("t", "trade")) {
  if (!file.exists(path)) {
    stop("Options CSV not found: ", path)
  }

  message("Reading: ", path)
  dt <- data.table::fread(
    path,
    select = keep_cols,
    nrows = if (is.finite(max_rows)) max_rows else -1
  )

  # Some exports (e.g. Spark) escape inner quotes as \" inside a quoted CSV field.
  # data.table::fread can often handle this automatically, but we defensively
  # unescape \" -> " so the JSON keys match our regex extraction.
  if ("trade" %in% names(dt)) {
    dt[, trade := gsub('\\"', '"', trade, fixed = TRUE)]
  }

  # timestamp in ms (prefer column t)
  dt[, timestamp_ms := as.numeric(t)]
  dt[, datetime_utc := as.POSIXct(timestamp_ms / 1000, origin = "1970-01-01", tz = "UTC")]
  dt[, date := as.Date(datetime_utc, tz = "UTC")]

  # Fast regex extraction from the JSON string in column 'trade'
  # We deliberately avoid jsonlite::fromJSON row-by-row (too slow for large files).
  tr <- dt[["trade"]]

  # instrument_name
  instr <- stringi::stri_match_first_regex(tr, '"instrument_name":"([^"]+)"')[, 2]
  iv    <- stringi::stri_match_first_regex(tr, '"iv":([0-9.]+)')[, 2]
  amt   <- stringi::stri_match_first_regex(tr, '"amount":([0-9.]+)')[, 2]
  idxp  <- stringi::stri_match_first_regex(tr, '"index_price":([0-9.]+)')[, 2]

  dt[, instrument_name := instr]
  dt[, iv := as.numeric(iv)]
  dt[, amount := as.numeric(amt)]
  dt[, index_price := as.numeric(idxp)]

  # parse instrument
  parsed <- parse_instrument_name(dt$instrument_name)
  dt <- cbind(dt, parsed)

  # basic cleaning
  dt <- dt[is.finite(iv) & is.finite(strike) & !is.na(expiry_date)]
  dt[, iv_dec := iv / 100]  # Deribit-style IV is typically in percent
  dt <- dt[iv_dec > 0 & iv_dec < 5]
  dt[, amount := fifelse(is.finite(amount) & amount > 0, amount, 1.0)]

  # keep only necessary columns
  dt <- dt[, .(
    date,
    datetime_utc,
    timestamp_ms,
    instrument_name,
    expiry_date,
    strike,
    cp_flag,
    iv_dec,
    amount,
    index_price
  )]
  data.table::setorder(dt, date)
  return(dt)
}
