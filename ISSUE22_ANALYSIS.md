# Issue #22: Statistical Analysis and Solutions

## Problem Statement

The `staggered` package fails with "non-conformable arguments" error when:
1. **Single cohort periods**: Only one treatment cohort has data at certain event times
2. **Out-of-range periods**: Event times fall outside available data range

## Statistical/Econometric Analysis

### Why Single-Cohort Periods are Problematic

#### Identification Requirement
DID estimation requires **at least two cohorts** at each event time to:
1. **Construct treatment-control comparisons**: Need a "treated" cohort and a "control" (not-yet-treated or never-treated) cohort
2. **Estimate variance**: Single observation (cohort) cannot reliably estimate sampling variability

#### Mathematical Explanation
The efficient estimator computes:
```
θ_hat = Σ_g A_θ,g * Y_bar_g - β' * Σ_g A_0,g * Y_bar_g
```

When only ONE cohort g has data at time t:
- Cannot form valid A_θ matrix (needs multiple cohorts for weighting)
- Matrix operations A[,gindices] fail when gindices has length 1
- No variance can be computed (need >= 2 observations)

### Scientifically Valid Solutions

#### ✅ Solution 1: Graceful Skip (RECOMMENDED)
**Principle**: "Only report what can be validly identified"

**Implementation**:
- When `length(eventTime) > 1`, try each event time separately
- Skip periods with insufficient cohort variation
- Return results only for identifiable periods
- Warn users about skipped periods

**Economic Validity**: ✓ 
- Does not force invalid estimates
- Transparent about limitations
- Common in empirical work (e.g., "We report estimates for periods with sufficient variation")

#### ✅ Solution 2: Use Alternative Control Groups
**Principle**: "Expand the set of valid comparisons"

**Options**:
a) **Use never-treated as controls** (`use_last_treated_only=FALSE`, default)
   - May help some periods gain additional cohorts
   
b) **Use last-treated-only** (`use_last_treated_only=TRUE`)
   - Following Sun & Abraham (2021)
   - May reduce identifiable periods but cleaner identification

**Economic Validity**: ✓
- Theoretically grounded (CS 2021, SA 2021)
- Changes the estimand slightly but remains interpretable

#### ⚠️ Solution 3: Pool Adjacent Periods
**Principle**: "Aggregate across time to gain sample size"

**Method**: Combine eventTime=-4,-3 into a single "pre-period" estimate

**Economic Validity**: ⚠️ Conditional
- Requires assumption: treatment effects constant across pooled periods
- Often too strong for event studies (we want period-specific effects!)
- Not recommended for `staggered` package

#### ❌ Solution 4: Imputation/Extrapolation
**Principle**: "Fill in missing cohort-time cells"

**Economic Validity**: ❌ NOT RECOMMENDED
- Introduces strong untestable assumptions
- Violates data-driven philosophy of modern DID
- Not done in credible empirical work

## Recommended Implementation for issue #22

### Primary Fix (Current PR)
1. **Handle single-cohort sapply issue**: Ensure sapply result is always a matrix
   - Fixes the "non-conformable arguments" error
   - Allows single-cohort periods to be estimated when possible

2. **Clear error messages**: When a period truly cannot be identified (out-of-range)
   - Informative error explaining why
   - Suggest valid eventTime range

### Future Enhancement (Separate PR)
Implement "graceful skip" for vector eventTime:

```r
# Pseudo-code for graceful skip
if (length(eventTime) > 1) {
  valid_results <- list()
  skipped_times <- c()
  
  for (et in eventTime) {
    result <- tryCatch({
      # Try to estimate for this event time
      staggered(df, estimand="eventstudy", eventTime=et, ...)
    }, error = function(e) {
      skipped_times <- c(skipped_times, et)
      NULL  # Skip this one
    })
    
    if (!is.null(result)) {
      valid_results[[length(valid_results) + 1]] <- result
    }
  }
  
  if (length(skipped_times) > 0) {
    warning(paste0("Skipped eventTime periods ", 
                   paste(skipped_times, collapse=", "),
                   " due to insufficient cohort variation"))
  }
  
  return(do.call(rbind, valid_results))
}
```

## Literature Support

### Callaway & Sant'Anna (2021)
- Emphasize identification requires "not-yet-treated" cohorts
- No valid comparison → no valid estimate

### Roth & Sant'Anna (2023) 
- Efficient estimator requires sufficient variation for β estimation
- "Asymptotically unbiased" assumes large sample (many cohorts)

### Sun & Abraham (2021)
- Recommend using "last-treated" to avoid contamination
- Acknowledge trade-off: cleaner identification vs. fewer identifiable periods

## Conclusion

**The single-cohort problem is NOT a bug** - it's a fundamental **identification limitation**.

**Correct approach**:
1. ✓ Fix sapply dimension issue (current PR)
2. ✓ Provide clear error messages
3. ✓ Document the limitation
4. ✓ (Future) Implement graceful skip for multiple eventTimes

**Do NOT**:
- ❌ Force estimates when identification fails
- ❌ Use imputation/extrapolation
- ❌ Hide the problem from users
