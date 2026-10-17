# Allele-aware drug repurposing for BUB1B mosaic variegated aneuploidy (MVA1)

Submission to **Rare Disease, Real Kid: MVA Hackathon 2026 - Track 2 (Drug Repurposing)**
(Sage Bionetworks, MVA Society, Hugging Face, BEACON). HF participant: `zbignevich1`.

> **Research hypotheses requiring clinical evaluation.** Nothing in this repository is medical advice or a treatment
> recommendation. It contains no doses or instructions for use in children.
>
> **No controlled hackathon data is included.** Submission 2 instantiates the framework for the proband: the provided VCF
> was processed locally on a separate machine by `src/vcf_panel_extract.py` (16-gene panel only) and then deleted.
> The repository holds only derived findings permitted by the organisers' data-use guidance: the two named *BUB1B*
> variants and their public annotations (`config/patient_genotype.yaml`) and two HPO terms used for drug safety.

## Idea in one paragraph

MVA1 alleles fall into two functional classes: truncating alleles make no transcript, while missense alleles in or
near the BUBR1 kinase domain make protein that is degraded fast; other missense alleles keep normal abundance but impair
checkpoint function (Suijkerbuijk et al., 2010). A drug that raises BUBR1 dosage can only help if one allele still makes
protein. The pipeline therefore classifies alleles, scores ten approved medicines on six referenced evidence
dimensions, ranks them separately for five biallelic genotype configurations, tests robustness over 5,000 perturbed
weightings and generates a do-no-harm list of checkpoint-dependent approved drugs.

Report: [`report/zbignevich1_track2_report.md`](report/zbignevich1_track2_report.md) (PDF alongside).

## Patient-specific result (Submission 2)

Proband genotype: *BUB1B* c.2210T>G p.(Leu737Ter) (ClinVar pathogenic nonsense, NMD-predicted) and c.3006T>G
p.(Asn1002Lys) (novel, pseudokinase domain, 10 residues from the unstable MVA1 allele L1012P), both heterozygous, phase
assumed trans. ESM-2, ThermoMPNN (AlphaFold model) and AlphaMissense place N1002K with destabilising alleles; ESM-1v
does not. Ranking with P(unstable) = 0.8 and HPO risk modifiers: niacin (rank 1 in 72% of 5,000 perturbed draws),
everolimus, dasatinib, metformin, sirolimus, arimoclomol (rank 2 if instability is confirmed). Decisive experiment:
BUBR1-N1002K half-life by cycloheximide chase. Details: report Section 2.

## Reproduce

```bash
pip install -r requirements.txt
python run_all.py                     # verify refs -> ClinVar landscape -> drug annotation -> scoring -> figures -> tests
python report/build_report.py --pdf   # rebuild the report (PDF rendering needs Edge, Chrome or Chromium)
```

All HTTP responses (PubMed, ClinVar, UniProt, ChEMBL, openFDA, Open Targets) are cached in `data/raw/http_cache`,
so reruns are deterministic and work offline.

Instantiate with a genotype (runs locally, sends nothing anywhere):

```bash
cd src && python classify_variants.py "c.2210T>G (p.Leu737Ter)" "c.3006T>G (p.Asn1002Lys)"
```

Patient layer (Submission 2):

```bash
python src/vcf_panel_extract.py --vcf /secure/WGS.vcf.gz --workdir /secure/work   # only where the controlled VCF lives
python src/variant_effect_models.py esm --workdir /tmp/w                          # public inputs; ~30 CPU-min
python src/variant_effect_models.py ddg --workdir /tmp/w --thermompnn /path/ThermoMPNN
python src/score.py && python src/patient_figures.py
```

## Layout

| Path | Content |
|---|---|
| `config/evidence.yaml` | Knowledge base: mechanism edges, allele classes, strategies, candidates, scores, do-no-harm patterns; every item carries citation keys |
| `config/references.yaml` | 68 references (PMID, DOI, title), verified against PubMed on every build |
| `config/weights.yaml` | Evidence weights, genotype configurations, sensitivity-analysis settings |
| `src/verify_refs.py` | PubMed verification; a missing or unverifiable citation stops the build |
| `src/clinvar_landscape.py` | Allele-class landscape for 7 MVA genes; kinase-window rule check |
| `src/drug_annotations.py` | ChEMBL / openFDA annotation, liability scan, Open Targets genes |
| `src/score.py` | Allele-aware scoring and Monte Carlo rank stability |
| `src/make_figures.py` | Figures in `figures/` |
| `src/classify_variants.py` | Genotype -> allele classes -> ranked candidates |
| `src/vcf_panel_extract.py` | Panel extraction + VEP annotation from a controlled VCF (run where the data live) |
| `src/variant_effect_models.py` | ESM-1v / ESM-2 scores, AlphaFold context, ThermoMPNN ddG (public inputs) |
| `src/patient_figures.py` | Patient figures (alleles, variant effect, class sweep) |
| `config/patient_genotype.yaml` | The proband's two named variants and evidence (derived findings only) |
| `logs/multimodel_critique/` | Verbatim log of the external-model critique round (no patient data in prompts) |
| `outputs/` | All result tables (CSV) |
| `tests/` | pytest: classification, citation completeness, genotype fit, determinism |
| `report/` | Report source, builder, Markdown and PDF |

## AI disclosure

AI use: Anthropic Claude Opus 5.5 via Claude Code, consumer Claude Max subscription, training setting OFF since
17 October 2026 (before any hackathon data were opened), no feedback or ratings given. Raw controlled files were
processed only by local scripts; Claude saw derived summaries (four panel variant calls with annotations, HPO terms).
VEP REST annotated 102 panel exon records. A critique round used free tiers of Gemini 3.7 Flash and Nemotron 3 Ultra
(four other calls failed) with prompts containing no patient data. Every reference is verified against PubMed (68/68).
Controlled data were deleted on 17 October 2026.

## License

Code: MIT (see `LICENSE`). Report, figures and result tables: CC BY 4.0, as required by the hackathon rules.

## Acknowledgement

This work was made possible through the Hackathon, organized by Sage Bionetworks in partnership with the MVA Society,
Hugging Face, and BEACON (The Benchmarking, Evaluation, and Assessment Consortium for Science), with prize sponsorship
from AWS and Anthropic. We are deeply grateful to the child and their family who generously contributed their data and
their story to advance research into this rare disease. We acknowledge their trust in making this Hackathon possible.
