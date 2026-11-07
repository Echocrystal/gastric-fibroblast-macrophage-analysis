# Figure QA: configuration overview

Scientific message: response distributions coexist with uneven MSI, EBV-ISH and CPS distributions across the four configurations.
Archetype: quantitative grid; A=response, B=MSI, C=EBV, D=CPS. Each panel has a distinct endpoint.
Author requested candidates directly; no additional figure-plan approval was needed.

- All 45 unique patient records included; fixed original group order and identity retained.
- Responses=12, MSI-positive=6, EBV-ISH-positive=5, CPS observed=41, CPS missing=4.
- MSI and EBV use the measured pentaplex/ISH assays, not WES estimates.
- No tests or performance metrics calculated. Wilson 95% intervals are descriptive binomial intervals.
- All 41 available CPS values plotted. Four missing values explicitly labelled; no imputation or negative classification.
- All bars originate at zero. Response and binary-composition ranges include 0–100%; CPS remains raw 0–100.
- Typography baseline copied; font floor 6.5 pt; Arial; top/right spines removed; legends outside data.
- Semantic status palette and existing Figure 6 phenotype colours, plus direct labels for colour redundancy.
- PDF/SVG preserve vector content and fonts; PNG preview at 300 dpi; exact canvas 183 × 165 mm.
- Asset scripts and companion PNGs inspected. They used incompatible paired/model structures, so their visual primitives were inherited without simulated data or inappropriate statistical overlays.
- Visual render inspected after layout repair: headings and notes separated; no data clipping; direct count labels legible. Repeated CPS zero values retain all individual points with deterministic vertical offsets, although dense zero stacks partially overlap.
- PDF page dimensions and PNG resolution checked: 183 × 165 mm, 300 dpi. No raw data values are changed by plotting offsets.
