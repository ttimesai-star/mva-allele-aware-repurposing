# Allele-aware drug repurposing for BUB1B mosaic variegated aneuploidy (MVA1)

Submission to **Rare Disease, Real Kid: MVA Hackathon 2026 - Track 2 (Drug Repurposing)**
(Sage Bionetworks, MVA Society, Hugging Face, BEACON). HF participant: `zbignevich1`.

> **Research hypotheses requiring clinical evaluation.** Nothing in this repository is medical advice or a treatment
> recommendation. It contains no doses or instructions for use in children.
>
> **No controlled hackathon data is included or was used.** Everything is built from public databases and literature.
> The framework is parametrised by genotype so the confirmed variants can be plugged in locally.

## Idea in one paragraph

MVA1 alleles fall into two functional classes: truncating alleles make no transcript, while missense alleles in or
near the BUBR1 kinase domain make protein that is degraded fast; other missense alleles keep normal abundance but impair
checkpoint function (Suijkerbuijk et al., 2010). A drug that raises BUBR1 dosage can only help if one allele still makes
protein. The pipeline therefore classifies alleles, scores ten approved medicines on six referenced evidence
dimensions, ranks them separately for five biallelic genotype configurations, tests robustness over 5,000 perturbed
weightings and generates a do-no-harm list of checkpoint-dependent approved drugs.

Report: [`report/zbignevich1_track2_report.md`](report/zbignevich1_track2_report.md) (PDF alongside).

## Reproduce

```bash
pip install -r requirements.txt
python run_all.py                     # verify refs -> ClinVar landscape -> drug annotation -> scoring -> figures -> tests
python report/build_report.py --pdf   # rebuild the report (PDF rendering needs Edge, Chrome or Chromium)
```

All HTTP responses (PubMed, ClinVar, UniProt, ChEMBL, openFDA, Open Targets) are cached in `data/raw/http_cache`,
so reruns are deterministic and work offline.

Instantiate with a confirmed genotype (runs locally, sends nothing anywhere):

```bash
cd src && python classify_variants.py "p.Gln383Ter" "p.Ile909Thr"
```

## Layout

| Path | Content |
|---|---|
| `config/evidence.yaml` | Knowledge base: mechanism edges, allele classes, strategies, candidates, scores, do-no-harm patterns; every item carries citation keys |
| `config/references.yaml` | 47 references (PMID, DOI, title), verified against PubMed on every build |
| `config/weights.yaml` | Evidence weights, genotype configurations, sensitivity-analysis settings |
| `src/verify_refs.py` | PubMed verification; a missing or unverifiable citation stops the build |
| `src/clinvar_landscape.py` | Allele-class landscape for 7 MVA genes; kinase-window rule check |
| `src/drug_annotations.py` | ChEMBL / openFDA annotation, liability scan, Open Targets genes |
| `src/score.py` | Allele-aware scoring and Monte Carlo rank stability |
| `src/make_figures.py` | Figures in `figures/` |
| `src/classify_variants.py` | Genotype -> allele classes -> ranked candidates |
| `outputs/` | All result tables (CSV) |
| `tests/` | pytest: classification, citation completeness, genotype fit, determinism |
| `report/` | Report source, builder, Markdown and PDF |

## AI disclosure

AI use: Anthropic Claude Opus 5.5 via Claude Code on a consumer Claude Max subscription. The account's 'Help improve our AI models' (training) setting was switched OFF on 17 October 2026, before submission; no hackathon-provided (gated) data was ever downloaded or shared with any model - the analysis used only public sources (PubMed, ClinVar, ChEMBL, openFDA, Open Targets). Models were used for literature triage, code and drafting; every reference was verified programmatically against PubMed (47/47).

## License

Code: MIT (see `LICENSE`). Report, figures and result tables: CC BY 4.0, as required by the hackathon rules.

## Acknowledgement

This work was made possible through the Hackathon, organized by Sage Bionetworks in partnership with the MVA Society,
Hugging Face, and BEACON (The Benchmarking, Evaluation, and Assessment Consortium for Science), with prize sponsorship
from AWS and Anthropic. We are deeply grateful to the child and their family who generously contributed their data and
their story to advance research into this rare disease. We acknowledge their trust in making this Hackathon possible.
