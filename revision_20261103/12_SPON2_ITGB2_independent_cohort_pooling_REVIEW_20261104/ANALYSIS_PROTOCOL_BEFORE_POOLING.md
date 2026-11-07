# Additional exploratory pooled analysis: frozen specification

Written on 4 November 2026 before computing the new pooled estimates/P values.
The individual-cohort results were already known when this question was proposed;
this is not a preregistered or independent validation analysis.

Question: Across the two independent spatial cohorts, is anchor-spot SPON2
positively associated with the exact ring-1 endpoint? Tests remain two-sided.

- Inputs: frozen section and donor coefficients in the independent 01 review
  package, generated with the original covariates, eligibility gates and fixed
  signatures. No raw spot-model refitting, expression imputation or smoothing.
- Three endpoints, all reported: ITGB2-plus-receiver composite, receiver score
  alone, ITGB2 alone. Receiver is the fixed primary C3-related signature.
- Included: GSE270678 main ST,19 donors; GSE251950 primary,9 patients.
  ST_R overlaps main donors and is excluded from pooling. All 28 individuals
  have all three endpoints. No selection using new P values.
- Verify frozen donor estimates against equal-section Fisher-z means from the
  accepted original sections. Verify original cohort coefficients/P/q values.
- Main weighting for THIS NEW analysis only: equal individuals; average the 28
  donor/patient Fisher-z coefficients and transform with tanh. Cohort weights
  are therefore19/28 and9/28. Do not weight by spot count.
- Prespecified sensitivity: each cohort50%; each individual within a cohort
  shares its cohort weight. This is not a random-effects meta-analysis.
- Pointwise95%CI: percentile bootstrap with99,999 replicates, resample19 and9
  individuals within their respective cohorts. Same resampling indexes for
  the three endpoints and both weighting schemes. Quantiles use linear
  interpolation. Fixed seed2026110412; no seed search.
- P: exhaustive two-sided individual-level sign-flip of Fisher-z coefficients,
  with absolute weighted sum at least as extreme as observed (tolerance1e-15).
  Count all2^28 assignments using a meet-in-the-middle algorithm. No Monte
  Carlo add-one correction. Assumes independent zero-centred symmetric
  individual coefficients under the sign-flip null; not assumption-free
  inference for an arbitrary cancelling average of heterogeneous effects.
- Multiplicity: BH across the3 pooled endpoint tests for the donor-equal
  scheme; separate3-test BH family for equal-cohort sensitivity. These new
  q values do not replace the original4-axis/3-supporting-endpoint families
  or correct the earlier candidate selection or the decision to pool.
- Algorithm validation: compare exact tail counts with brute force for small
  synthetic vectors and reproduce original per-cohort sign-flip P values.
- Descriptive sensitivity: report leave-one-individual-out pooled effect
  range with weights recalculated under each scheme; no new LOPO P values.
- Report all original cohort results, both pooling schemes and full individual
  source data regardless of significance. No claim of component superiority,
  fractional contribution, mediation, direct binding or causal signalling.
- Deliverables: independent scripts, TSV/JSON/XLSX tables, review forest plot,
  standalone Word/PDF report, input hashes and verification record.
- Existing figures, Word documents, formal supplemental tables and GitHub
  archive remain unchanged. This folder alone receives new outputs.
