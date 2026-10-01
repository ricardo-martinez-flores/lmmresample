# lmmresample — interface specification (v0.1)

Status: draft for review. This document defines the public interface before
implementation. Function names, arguments and return values here are binding
for v0.1 once approved; internal implementation details are not.

---

## 1. Scope

### 1.1 What the package does

Resampling-based inference (permutation and bootstrap) for a single model
fitted to repeated measurements organised in **blocks**: trials, epochs, task
blocks or participants. The resampling unit is always a complete block, never
an individual sample. Each replicate refits the user's model on the full data,
so any dependence within a block (e.g. autocorrelation across samples of a
trial) is carried intact into the reference distribution without having to be
modelled explicitly.

### 1.2 Signals

Any densely sampled signal with a block structure, including:

- oculomotor: pupil size, saccades, gaze position;
- neurophysiological: EEG/ERP, MEG, fNIRS (HbO/HbR epochs);
- autonomic: heart rate and HRV, electrodermal activity, respiration;
- motor and biomechanical: EMG, kinematics, force plate, accelerometry;
- continuous behavioural measures: mouse or joystick tracking, continuous
  ratings.

### 1.3 Levels of aggregation

The same functions apply at every level; only the fitted model and the
declared exchangeability change.

| Level | Typical model | Permutation unit | Bootstrap unit |
|---|---|---|---|
| Full time series (samples within trials) | `lmer`, many rows per trial | trial, within participant | participant |
| Feature per trial (peak, latency, area) | `lmer`, one row per trial | trial, within participant | participant |
| Feature per participant × condition | `lmer`, random intercept | condition label, within participant | participant |
| Feature per participant | `lm` / `glm` | group label, between participants | participant |

### 1.4 Supported models

v0.1: `stats::lm`, `stats::glm`, `lme4::lmer`, `lme4::glmer`.
Later: `glmmTMB::glmmTMB`.

Refitting uses `stats::update()` with a modified `data` argument, so the
package relies only on the model object exposing its call and data.

### 1.5 Requirements and limits (to be stated in the documentation)

- **Block structure is required.** A single continuous recording per
  participant with no trial structure (e.g. a 24-h Holter) supports only
  between-participant comparisons, with the participant as the block.
  Resampling within one long series (temporal block bootstrap) is out of scope.
- **Carry-over between consecutive trials** (e.g. slow haemodynamic responses
  with short inter-trial intervals) violates exchangeability of trials.
  `perm_calibrate()` is the recommended way to detect this.
- The tested variable must be constant within the permutation unit (e.g.
  condition is constant within a trial).

---

## 2. Design principles

1. **Model first.** The user fits the model as usual; every function takes the
   fitted model as its first argument.
2. **Explicit exchangeability.** The user always declares what is exchanged.
   There is no silent default that could be wrong for the design.
3. **Store everything.** Result objects keep the full null or bootstrap
   distribution, refit failures, the seed and the exchangeability
   specification, so results can be inspected, re-summarised and reproduced.
4. **Consistent arguments.** The same argument names mean the same thing in
   every function (`term`, `unit`, `exchange`, `cluster`, `B`, `seed`).
5. **Parallel by configuration, not by argument.** Computation uses the
   `future` framework; the user chooses the backend with `future::plan()`.
   Results are reproducible regardless of the backend.
6. **Fail loudly, degrade gracefully.** Invalid designs stop with a clear
   message. Individual refits that fail to converge are recorded and reported,
   not silently dropped.

---

## 3. Function overview

Functions are grouped by prefix so that they are easy to find with
autocompletion.

| Family | Function | Purpose |
|---|---|---|
| Permutation | `perm_test()` | Permutation test for one term |
| | `perm_maxt()` | Family-wise error control across several terms or models (max-t) |
| | `perm_calibrate()` | Type I error calibration of the permutation test by simulation |
| Bootstrap | `boot_ci()` | Bootstrap confidence intervals for fixed effects |
| Diagnostics | `diag_loso()` | Leave-one-subject-out influence |
| | `diag_agreement()` | Model coefficient vs. mean participant-level difference |
| | `diag_acf()` | Within-series autocorrelation of residuals |
| Exchangeability | `exch_within()`, `exch_free()`, `exch_custom()` | Declare what is exchanged |
| Simulation | `sim_blocks()` | Simulate block-structured data for examples, power and calibration |

Note: earlier drafts used `perm_lmer()`, `boot_lmer()`, `maxt_lmer()` and
`loso_lmer()`. The names were generalised because `lm` and `glm` models are
now supported.

---

## 4. Exchangeability

```r
exch_within(block)            # permute units within each level of `block`
exch_free(strata = NULL)      # permute units freely, optionally within strata
exch_custom(how)              # any permute::how() specification (advanced)
```

Together with `unit` (the column identifying the resampled block), these
cover the standard designs:

| Design | Tested variable | `unit` | `exchange` |
|---|---|---|---|
| Within-participant | condition | `"trial"` | `exch_within("participant")` |
| Between-group | group | `"participant"` | `exch_free()` |
| Between-group, stratified | group | `"participant"` | `exch_free(strata = "site")` |
| Repeated measures (participant × condition rows) | condition | `"row"` or condition-level id | `exch_within("participant")` |

Validation performed before any refit:

- `unit` and the block/strata columns exist in the data;
- the tested variable is constant within each `unit`;
- for `exch_within()`, each block contains more than one level of the tested
  variable (otherwise permutation is impossible in that block, and the user is
  warned how many blocks are uninformative);
- the number of distinct permutations is reported; if it is smaller than `B`,
  the user is warned and an exact test is suggested.

---

## 5. Function reference

### 5.1 `perm_test()`

```r
perm_test(
  model,
  term,
  unit,
  exchange,
  coef = NULL,
  alternative = c("two.sided", "greater", "less"),
  B = 1999,
  seed = NULL,
  data = NULL
)
```

- `term`: name of the variable to permute (a column of the data).
- `coef`: name of the model coefficient used as the statistic. Inferred when
  `term` maps to a single coefficient (numeric or two-level factor); required
  otherwise. Multi-level factors (F / likelihood-ratio statistics) are planned
  for a later version.
- Statistic: the t (or z for `glm`/`glmer`) value of `coef`.
- p-value: Phipson–Smyth, *p* = (*b* + 1) / (*B*′ + 1), where *B*′ is the
  number of successful refits.
- `data`: defaults to the data stored in the model call; needed only if the
  data cannot be recovered from the model.

Returns an object of class `lmmr_perm` containing: observed statistic,
null distribution, p-value, `B` requested and successful, number of
non-convergent and singular refits, seed, exchangeability specification and
call.

Methods: `print()`, `summary()`, `tidy()`, `plot()`.

### 5.2 `perm_maxt()`

```r
perm_maxt(
  ...,                        # named perm_spec() objects
  unit,
  exchange,
  method = c("single-step", "step-down"),
  B = 1999,
  seed = NULL
)

perm_spec(model, term, coef = NULL)
```

- Each `perm_spec()` pairs a model with the term tested in it. Several terms
  from the same model or the same term across several models (e.g. one model
  per feature) are both allowed.
- All specifications must share the permutation units: the **same shuffle is
  applied to every model in each replicate**, which preserves the dependence
  among statistics.
- Statistic: |t| (two-sided). Adjusted p-values by the Westfall–Young max-t
  procedure, single-step or step-down.

Returns `lmmr_maxt`: observed |t| per specification, unadjusted and adjusted
p-values, the joint null matrix (B × number of specifications) and refit
diagnostics per model.

### 5.3 `perm_calibrate()`

```r
perm_calibrate(
  model,
  term,
  unit,
  exchange,
  coef = NULL,
  null = c("sharp", "mean"),
  ar1 = NULL,
  series = NULL,
  time = NULL,
  n_sim = 200,
  B = 199,
  alpha = 0.05,
  seed = NULL
)
```

Simulates data under H0 from the user's own **reduced model** (the model with
`term` removed), runs `perm_test()` on each simulated dataset and reports the
distribution of p-values.

- `null = "sharp"`: the effect of `term` is zero for every participant (random
  slopes for `term` are also removed). This is the null hypothesis that a
  permutation test is exact for.
- `null = "mean"`: the average effect is zero but participants vary around it
  (random slope variance kept). Shows how the test behaves under the weaker,
  often more realistic null.
- `ar1`: optional AR(1) coefficient for the residuals within each `series`
  ordered by `time`. Can be a number or the output of `diag_acf()`.

Returns `lmmr_calib`: simulated p-values, rejection rate at `alpha` with a
binomial confidence interval, and a uniformity test.

Cost: `n_sim × B` refits. Defaults are chosen for a quick check; the
documentation will give guidance for a definitive one.

### 5.4 `boot_ci()`

```r
boot_ci(
  model,
  cluster,
  resample = c("case", "wild"),
  ci = c("percentile", "bca", "basic"),
  level = 0.95,
  terms = NULL,
  B = 1999,
  seed = NULL,
  data = NULL
)

confint(object, type = c("percentile", "bca", "basic"), level = 0.95)
```

- `resample = "case"`: resample whole clusters (participants) with
  replacement; resampled copies receive new identifiers so the model treats
  them as distinct participants. Nested designs only in v0.1.
- `resample = "wild"`: residual-based wild bootstrap; robust to
  heteroscedasticity and suitable for few clusters.
- `ci`: interval type. All replicates are stored, so `confint()` can return
  another interval type afterwards without refitting.
- BCa acceleration is estimated by jackknife over clusters (leave one cluster
  out). These refits are the same ones used by `diag_loso()` and are reused
  when available.
- `terms`: fixed effects to report; all by default.
- Case resampling requires nested random effects. If the model contains
  crossed random effects (e.g. participants and items), `boot_ci()` warns that
  the intervals reflect sampling of clusters only, and points to
  `lme4::bootMer()` as an alternative.

Returns `lmmr_boot`.

### 5.5 `diag_loso()`

```r
diag_loso(model, cluster, terms = NULL, data = NULL)
```

Refits the model leaving out each cluster in turn. Returns the coefficients
for each omission and their change relative to the full model.

### 5.6 `diag_agreement()`

```r
diag_agreement(model, term, cluster, tolerance = 0.25, data = NULL)
```

Compares the model coefficient for a two-level within-cluster `term` with the
mean of the participant-level differences computed directly from the data.
Warns when the relative discrepancy exceeds `tolerance`. A frequent cause in
time-series models is a missing fixed effect of time when time courses differ
between conditions; the message suggests checking this.

### 5.7 `diag_acf()`

```r
diag_acf(model, series, time, lag_max = NULL)
```

Autocorrelation of the model residuals within each `series` (e.g. trial),
ordered by `time`, averaged across series. Returns the ACF and an AR(1)
estimate that can be passed to `perm_calibrate(ar1 = )`.

### 5.8 `sim_blocks()`

```r
sim_blocks(
  design = c("within", "between", "mixed"),
  level = c("timeseries", "trial", "participant"),
  n_participants = 30,
  n_trials = 40,
  n_time = 100,
  effect = 0,
  ar1 = 0.9,
  ...
)
```

Generates block-structured data with known effects. Used in examples, tests
and vignettes, and useful to users for power analysis.

---

## 6. Plots

Every result object has a `plot()` method returning a `ggplot` object, so
users can modify it. `ggplot2` is a suggested dependency; plotting functions
stop with an informative message if it is not installed. The data behind each
plot is available through `tidy()`.

| Object | `plot(type = )` | Shows |
|---|---|---|
| `lmmr_perm` | `"null"` (default) | Null distribution with the observed statistic |
| | `"trace"` | Running Monte Carlo p-value as B increases |
| `lmmr_maxt` | `"null"` (default) | Null distribution of max\|t\| with each observed \|t\| and the critical value |
| `lmmr_calib` | `"pp"` (default) | Simulated p-values against the uniform, with confidence band |
| | `"rejection"` | Rejection rate at α with its interval |
| `lmmr_boot` | `"dist"` (default) | Bootstrap distribution per term with estimate and intervals (percentile and BCa overlaid) |
| `lmmr_loso` | `"influence"` (default) | Coefficient when each cluster is left out, with the full estimate and its interval |
| `lmmr_agreement` | `"clusters"` (default) | Participant-level differences against the model coefficient |
| `lmmr_acf` | `"acf"` (default) | Mean residual ACF with spread across series |

Generic residual diagnostics (residuals vs fitted, QQ) are left to existing
packages such as `performance` and `DHARMa`.

---

## 7. Common behaviour

- **Reproducibility.** `seed` sets the RNG for the whole call, including
  parallel workers (`future.seed`). The same seed gives the same result with
  any `future` plan.
- **Progress.** Reported through `progressr`, which the user can enable or
  silence.
- **Refit failures.** Non-convergent and singular refits are counted per
  call. By default non-convergent refits are excluded and the effective `B` is
  reported; `print()` warns when more than 5% of refits failed.
- **Messages.** Errors and warnings use `cli`, naming the argument at fault
  and how to fix it.

---

## 8. Usage examples

```r
library(lmmresample)
library(lme4)
```

### 8.1 Full time series, within-participant condition

```r
d <- sim_blocks("within", level = "timeseries", effect = 0.1)
m <- lmer(pupil ~ condition + time + (1 + condition | participant), data = d)

diag_acf(m, series = "trial", time = "time")
diag_agreement(m, term = "condition", cluster = "participant")

cal <- perm_calibrate(m, term = "condition", unit = "trial",
                      exchange = exch_within("participant"),
                      ar1 = diag_acf(m, "trial", "time"),
                      series = "trial", time = "time")
plot(cal)

p <- perm_test(m, term = "condition", unit = "trial",
               exchange = exch_within("participant"))
b <- boot_ci(m, cluster = "participant", ci = "bca")
plot(p); plot(b)
```

### 8.2 Feature per trial, several features with family-wise control

```r
d <- sim_blocks("within", level = "trial")
m_peak    <- lmer(peak    ~ condition + (1 | participant), data = d)
m_latency <- lmer(latency ~ condition + (1 | participant), data = d)

perm_maxt(
  peak    = perm_spec(m_peak,    "condition"),
  latency = perm_spec(m_latency, "condition"),
  unit = "trial", exchange = exch_within("participant")
)
```

### 8.3 Between-group comparison, one value per participant

```r
d <- sim_blocks("between", level = "participant")
m <- lm(feature ~ group + age, data = d)

perm_test(m, term = "group", unit = "participant", exchange = exch_free())
boot_ci(m, cluster = "participant", ci = "percentile")
```

### 8.4 Mixed design, between-group effect on the full time series

```r
d <- sim_blocks("mixed", level = "timeseries")
m <- lmer(pupil ~ group * condition + time + (1 + condition | participant),
          data = d)

perm_test(m, term = "group", unit = "participant", exchange = exch_free())
perm_test(m, term = "condition", unit = "trial",
          exchange = exch_within("participant"))
```

### 8.5 Parallel computation

```r
future::plan("multisession", workers = 8)
perm_test(m, term = "condition", unit = "trial",
          exchange = exch_within("participant"), B = 9999, seed = 1)
```

---

## 9. Dependencies

- Imports: `stats`, `lme4`, `future.apply`, `cli`, `rlang`.
- Suggests: `ggplot2`, `progressr`, `permute` (for `exch_custom()`),
  `lmeresampler` (comparisons in tests), `testthat`, `knitr`, `rmarkdown`.

---

## 10. Roadmap

Development milestones. The first CRAN release will be v0.2, which is also the
version described in the methods paper. v0.1 is an internal milestone, tagged
on GitHub only.

**v0.1 — core (internal milestone)**
Everything specified in sections 3–8.

**v0.2 — covariates, interactions and more models (first CRAN release; paper)**
- Freedman–Lane and ter Braak permutation with nuisance covariates, which also
  provides valid tests of interactions (e.g. `group × condition` in mixed
  designs). Calibration is mandatory in the documentation of these methods.
- F / likelihood-ratio statistics for multi-level factors and their
  interactions, in `perm_test()` and `perm_maxt()`.
- `glmmTMB` models.

**Not planned**
- Temporal block bootstrap within a single long series without block
  structure.

---

## 11. Open decisions

1. ~~Interactions in mixed designs.~~ Resolved: tested in v0.2 through
   Freedman–Lane. In v0.1, `perm_test()` stops with an informative message when
   `term` is involved in an interaction.
2. **Default `null` in `perm_calibrate()`.** `"sharp"` matches what the
   permutation test guarantees; `"mean"` is closer to typical data and more
   informative about real-world behaviour.
3. **Default `B`.** 1999 for tests and intervals (adequate for α = 0.05 and
   BCa) versus 4999 (more stable BCa endpoints at higher cost).
4. **Wild bootstrap implementation.** Build on `lmeresampler` (fewer lines,
   one more dependency) or implement directly (full control, BCa compatibility).
5. **Output class prefix.** `lmmr_*` is short; `lmmresample_*` is explicit.
