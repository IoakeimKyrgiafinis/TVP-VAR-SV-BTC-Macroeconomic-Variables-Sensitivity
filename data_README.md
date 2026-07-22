# Data

Raw data files are **not committed** to this repository. All series are fetched programmatically in `R/btc_tvpvar_full.R` at runtime.

## Sources

| Variable | Source | Ticker / Series ID | Frequency | Notes |
|----------|--------|--------------------|-----------|-------|
| Bitcoin (BTC-USD) | Yahoo Finance via `quantmod` | `BTC-USD` | Daily | Available from ~July 2014 |
| CBOE VIX | Yahoo Finance via `quantmod` | `^VIX` | Daily | |
| S&P 500 | Yahoo Finance via `quantmod` | `^GSPC` | Daily | |
| Gold ETF | Yahoo Finance via `quantmod` | `GLD` | Daily | SPDR Gold Shares |
| US Dollar Index | Yahoo Finance via `quantmod` | `DX-Y.NYB` | Daily | |
| EPU Index | FRED via `fredr` | `USEPUINDXD` | Daily | Baker, Bloom & Davis (2016) |

## FRED API Key

The EPU series is fetched from the Federal Reserve Economic Data (FRED) API. A free API key is required:

1. Register at [https://fred.stlouisfed.org/docs/api/api_key.html](https://fred.stlouisfed.org/docs/api/api_key.html)
2. Replace the placeholder key in line 6 of `R/btc_tvpvar_full.R`:
   ```r
   fredr_set_key("YOUR_API_KEY_HERE")
   ```

## Sample Period

- **Start:** 2014-07-18 (earliest reliable weekly BTC data)
- **End:** 2025-05-27
- **Frequency:** Weekly (last observation of each calendar week)
- **Usable observations after training sample and lag adjustment:** ~436 weeks

## Pre-processing

1. All daily series are downloaded and forward-filled (`na.locf`) to handle non-synchronous trading calendars.
2. Series are merged on a common daily index, then sampled at weekly frequency by taking the last observation of each week (`endpoints(..., on = "weeks")`).
3. Log-returns are computed: `diff(log(price))`.
4. All variables are standardised (zero mean, unit variance) prior to TVP-VAR estimation.
5. A training sample of τ = 120 weeks is used to initialise prior distributions and is excluded from the reported estimation window.
