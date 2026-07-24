# Bitcoin Market Integration and Macroeconomic Sensitivity
### A Time-Varying Parameter VAR Model with Stochastic Volatility

**Weekly Data: July 2014 – May 2025 · Six-Variable System: EPU · VIX · DXY · S&P 500 · Gold · Bitcoin**

---

## Overview

This project examines the time-varying relationship between Bitcoin (BTC) and five macro-financial variables using a six-variable **TVP-VAR-SV** (Time-Varying Parameter Vector Autoregression with Stochastic Volatility) model estimated on weekly data from July 2014 to May 2025.

The central question: has Bitcoin evolved from an isolated speculative asset into an integrated component of the global financial system — and if so, *when* and *through which channels*?

The model is identified via Cholesky decomposition of the time-varying covariance matrix. Inference is conducted through Bayesian MCMC (30,000 draws, 15,000 burn-in) using the [`bvarsv`](https://cran.r-project.org/package=bvarsv) package in R.

This paper is inspired by Panagiotidis, Stengos, and Vravosinos (2019), who examine the effects of stock market returns, exchange rates, gold, oil, central bank rates, and internet search intensity on bitcoin returns using standard VAR and Factor-Augmented VAR (FAVAR) models.This paper builds on their work by allowing the transmission mechanism to vary continuously over time via a TVP-VAR-SV specification, enabling us to trace how bitcoin's sensitivity to macroeconomic financial shocks has evolved across different market regimes, something that static VAR and FAVAR models are not designed to capture.
---

## Key Findings

| Driver | FEVD Share | Direction | Strongest Era | Notes |
|--------|-----------|-----------|---------------|-------|
| **VIX** | Up to 25% | Negative | COVID 2020 | Robust, persistent throughout |
| S&P 500 | Up to 8% | Positive → Negative | 2020–2022 | Regime change post-2022 |
| Gold | < 2.5% | Positive | COVID era | Stable but small |
| DXY | < 2.5% | Negative (h=1) | No clear peak | Weak, ambiguous |
| EPU | < 2.5% | Negative | 2017 Retail Era | Fades over time |

**VIX is the dominant external driver.** Bitcoin's self-shock still explains 60–100% of its forecast error variance, but financial market stress transmits consistently and negatively at the weekly structural frequency — an effect that deepened sharply during COVID and has not reverted.

### Three-Era Narrative

- **Retail Era (pre-2018):** Bitcoin was genuinely isolated. IRF responses to all external variables were statistically indistinguishable from zero. Dynamics dominated by idiosyncratic supply-demand cycles and exchange infrastructure shocks.
- **COVID Crisis Era (2020):** The March 2020 "Black Thursday" crash pulled Bitcoin into the global risk-off dynamic for the first time in a statistically robust sense. VIX FEVD contribution surged to ~25%; cumulative VIX response reached approximately −1.0.
- **Maturity Era (post-2022):** VIX sensitivity has moderated but remains elevated relative to 2017. The S&P 500 relationship underwent a complete sign reversal in the cumulative IRF, coinciding with the onset of Federal Reserve tightening and the FTX collapse.

---

## Data

All data is fetched programmatically — no raw files need to be committed. See `data/README.md` for full details.

| Variable | Description | Source | Ticker |
|----------|-------------|--------|--------|
| BTC | Bitcoin / USD | Yahoo Finance | `BTC-USD` |
| VIX | CBOE Volatility Index | Yahoo Finance | `^VIX` |
| S&P 500 | US equity benchmark | Yahoo Finance | `^GSPC` |
| Gold | SPDR Gold Shares ETF | Yahoo Finance | `GLD` |
| DXY | US Dollar Index | Yahoo Finance | `DX-Y.NYB` |
| EPU | Economic Policy Uncertainty Index | FRED | `USEPUINDXD` |

Raw daily series are forward-filled, merged to a common calendar, then sampled at weekly frequency (last observation of each week). All variables are transformed to log-returns and standardised prior to estimation.

> **FRED API key required.** Register for free at [fred.stlouisfed.org](https://fred.stlouisfed.org/docs/api/api_key.html) and replace the key in `R/btc_tvpvar_full.R` (line 6).

---

## Reproducing the Analysis

### 1. Install dependencies

```r
install.packages(c("quantmod", "fredr", "xts", "bvarsv", "vars", "coda"))
```





The script runs end-to-end in numbered steps:

| Step | Description |
|------|-------------|
| 1 | Download data from Yahoo Finance and FRED |
| 2–3 | Daily merge → weekly sampling → log-returns |
| 4 | Lag selection via information criteria (AIC/BIC/HQ) |
| 5 | Standardise variables |
| 6 | **Estimate TVP-VAR-SV** via Bayesian MCMC — saves `fit_tvp_6var.rds` |
| 7 | Dimension checks and NaN diagnostics |
| 8 | Geweke convergence diagnostics and inefficiency factors |
| 9 | Time-varying FEVD computation and plots (Figs 1–6) |
| 10 | Point-in-time IRFs with 68% / 95% credible bands (Figs 7–21) |
| 11 | Continuous IRF timelines with 68% bands (Figs 22–26) |
| 12 | Cumulative IRFs — point-in-time and continuous (Figs 27–46) |
| — | Export all plots to `BTC_TVP_VAR_Plots.pdf` |

> **Runtime note:** The MCMC estimation (Step 6) takes approximately **2–4 hours** on a modern laptop for 30,000 iterations.
> ```r
> fit_tvp <- readRDS("output/fit_tvp_6var.rds")
> ```
> 

---

## Model Specification

The TVP-VAR-SV model follows Primiceri (2005):

$$y_t = c_t + A_{1,t} y_{t-1} + A_{2,t} y_{t-2} + \varepsilon_t, \quad \varepsilon_t \sim \mathcal{N}(0, \Sigma_t)$$

Time-varying parameters evolve as random walks. The covariance matrix $\Sigma_t$ is decomposed via time-varying Cholesky factorisation with stochastic volatility. Lag order $p = 2$ selected by information criteria.

**Cholesky ordering (macro → financial → crypto):**

```
EPU → VIX → DXY → S&P 500 → Gold → BTC
```

This identifies Bitcoin as responding contemporaneously to all other shocks but not driving them contemporaneously — a standard assumption in the TVP-VAR literature for a small, fast-moving asset within a macro-financial system.

**MCMC settings:** 30,000 draws · 15,000 burn-in · training sample τ = 120 weeks (~2.3 years)

---

## Figures

All 46 figures are exported to `output/BTC_TVP_VAR_Plots.pdf`. Individual PNGs are in `output/figures/`.

| Figure range | Content |
|---|---|
| 1–6 | Time-varying FEVD: share of BTC variance explained by each shock |
| 7–21 | Point-in-time IRFs across three eras (Retail / COVID / Maturity) |
| 22–26 | Continuous IRF timelines (h=1) with 68% credible bands |
| 27–41 | Point-in-time cumulative IRFs (H=12) across three eras |
| 42–46 | Continuous cumulative IRF timelines (H=12) |

---

## References

- Primiceri, G. E. (2005). Time varying structural vector autoregressions and monetary policy. *Review of Economic Studies*, 72(3), 821–852.
- Krueger, F. (2015). bvarsv: Bayesian Analysis of a Vector Autoregressive Model with Stochastic Volatility and Time-Varying Parameters. R package v1.1.
- Panagiotidis, T., Stengos, T., & Vravosinos, O. (2019). The effects of markets, uncertainty and search intensity on bitcoin returns. *International Review of Financial Analysis*, 63, 220–242.
- Baker, S. R., Bloom, N., & Davis, S. J. (2016). Measuring economic policy uncertainty. *Quarterly Journal of Economics*, 131(4), 1593–1636.
- Geweke, J. (1992). Evaluating the accuracy of sampling-based approaches to the calculation of posterior moments. In *Bayesian Statistics 4*. Oxford: Clarendon Press.

---

## License

This project is released under the [MIT License](LICENSE). The working paper (`paper/`) is shared for academic reference.
