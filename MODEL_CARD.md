# Model Card: PVRIG-PortfolioRank v2.0

## Model summary

PVRIG-PortfolioRank is a deterministic, multi-objective computational pipeline
for generating and prioritizing VHH CDR designs intended to block the
PVRIG–PVRL2 interaction.

It is not a single end-to-end neural-network checkpoint. The released model
artifact is an auditable JSON configuration combining:

1. an IMGT-position-conditioned CDR residue prior;
2. PVRIG interface-composition constraints from 8X6B and 9E6Y;
3. an upstream ESM2 embedding model;
4. ANARCII2 antibody-likeness;
5. sequence-level developability filters;
6. diversity-aware selection and experimental portfolio constraints.

The exact task-specific configuration is
[`models/pvrig_portfolio_rank_v2.json`](models/pvrig_portfolio_rank_v2.json).
The upstream ESM2 weights are referenced by model ID and are not redistributed.

## Intended use

- retrospective reproduction of the released 50-candidate result;
- research comparison of CDR generation and prioritization strategies;
- prioritization of candidates for expression, binding and competition assays;
- software and scoring audit.

## Out-of-scope use

- clinical use;
- estimating a physical Kd, IC50 or blocking percentage from the released score;
- asserting therapeutic efficacy or safety;
- asserting freedom to operate or patentability;
- treating PVRIG–PVRL2 structures as candidate-antibody complex predictions.

## Inputs

- the 48-reference PVRIG positive library bundled with the organizer validator;
- PDB 8X6B and 9E6Y;
- five VHH frameworks contained in the positive library;
- the ESM2 model `facebook/esm2_t6_8M_UR50D`;
- random seed `20260726`.

## Stage 1: candidate generation and screening

### IMGT-aware CDR prior

ANARCII2 numbers all reference variable domains under IMGT. Each CDR is mapped
to 12 relative-position bins. Smoothed residue counts form a categorical prior.
Sampling introduces target-conditioned enrichment for polar/aromatic contacts
and modest acidic complementarity.

### Interface-composition prior

PVRIG residues within 5 Å of PVRL2 are extracted from 8X6B and 9E6Y. Generated
CDRs are scored by their deviation from a target composition over aromatic,
polar, acidic, basic, hydrophobic and Gly/Pro groups.

This is a composition prior only. The released run did not perform
candidate-specific antibody–PVRIG docking.

### ESM2 enrichment prior

ESM2 mean-pooled embeddings are calculated for positives, random decoys and
generated candidates. The direction from the decoy centroid to the positive
centroid defines a one-class enrichment axis. Candidate projections are
quantile-scaled within the screened batch.

### Hard filters

The CDR design region rejects:

- non-standard amino acids;
- added Cys or Met;
- N-X-S/T, NS, NG, DG, DS and DP motifs;
- four or more consecutive hydrophobic residues;
- triple residue repeats;
- excessive CDR charge, hydrophobicity or aromaticity;
- corresponding CDR identity at or above the design ceiling.

These filters do not imply that the full framework lacks all chemical
liabilities. Full-length framework motifs are reported separately.

### Stage-1 score

The first-pass score is:

```text
0.28 × ESM2 enrichment prior
+ 0.24 × interface-composition prior
+ 0.18 × CDR statistical prior
+ 0.12 × ANARCII2 antibody-likeness
+ 0.18 × sequence developability
```

The score is used only for within-batch prioritization.

## Stage 2: calibration and portfolio selection

The five stage-one components are converted to within-batch percentiles and
combined with a positive-reference novelty margin and an H3-length preference.
Two thousand Dirichlet weight perturbations produce a weight-sensitivity range.
This range is not an affinity confidence interval.

The first experimental portfolio uses framework quotas:

- F1: 2;
- F3: 4;
- F4: 2;
- F5: 2;
- at most one H3 longer than 18 residues.

The quota prevents the ten experimental slots from sharing one dominant
framework/loop regime.

## Validation

The released result was checked with the organizer's
`clickmab-bio/ab-data-validator` at commit
`97df17aa09bc576a861cf0d8242de97af379fd80`, using:

- ANARCII 2.0.8 / IMGT;
- global pairwise CDR alignment;
- identity threshold 0.75.

Result: 50/50 candidates passed; maximum observed corresponding-CDR identity
was 71.43%.

## Limitations and known risks

1. The positive set is small and patent-derived, so priors may over-represent
   known sequence families.
2. No candidate-specific complex geometry, explicit solvent binding free
   energy or molecular dynamics was used in the released score.
3. ESM2 enrichment is a one-class sequence prior, not a supervised affinity
   predictor.
4. Full-length expression, SEC purity, thermal stability, nonspecific binding
   and immunogenicity require experimental testing.
5. Accelerator-specific floating-point kernels may cause small full-regeneration
   ranking differences. The fixed Stage-1 snapshot supports exact replay.

## Required wet-lab follow-up

At minimum:

- expression yield and monomeric purity;
- SEC and DSF;
- direct PVRIG binding;
- PVRIG–PVRL2 competition/HTRF blocking;
- concentration-response confirmation;
- counterscreens against related nectin-axis receptors;
- forced-degradation testing for exposed Asn/Asp/Met liabilities.
