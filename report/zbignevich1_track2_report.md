# Allele-aware drug repurposing for BUB1B mosaic variegated aneuploidy (MVA1)

*Rare Disease, Real Kid: MVA Hackathon 2026 - Track 2 (Drug Repurposing) - submission 2 by zbignevich1, instantiated for the proband's BUB1B genotype*


> **Nature of this document.** This report proposes *research hypotheses requiring clinical evaluation*. It is not
> medical advice and not a treatment recommendation for any individual. It contains no doses, schedules or
> instructions for use in children. Every candidate below is an approved medicine proposed for laboratory and
> clinical evaluation by qualified investigators and the treating team; several carry serious toxicities.
> **What changed in this submission (2 of 3):** the controlled hackathon VCF and HPO profile were analysed locally,
> the framework was instantiated for the proband's confirmed-candidate genotype (Section 2), and the missense allele
> was assessed with sequence and structure models. Only the two named variants, their public annotations and the HPO
> terms relevant to drug safety are reported, as the organisers' data-use guidance allows; no genotype-scale data
> leave the analysis machine, and all copies were deleted after the analysis (Section 2.1).

**AI disclosure (required).** AI use: Anthropic Claude Opus 5.5 via Claude Code on a consumer Claude Max subscription, with the 'Help improve our AI models' (training) setting switched OFF on 17 October 2026, before any hackathon data were opened, and no feedback or ratings given on any output - the configuration the organisers' 28 September 2026 guidance accepts for consumer tiers. Raw controlled files (VCF, phenotype document) were processed only by local scripts on a separate workstation; Claude received only derived summaries (the four rare protein-altering calls in a 16-gene panel with zygosity and public annotations, per-gene counts, and the HPO table). Panel exon records (102) were annotated with Ensembl VEP REST, a processing tool. A critique round used free tiers of Google Gemini 3.7 Flash and NVIDIA Nemotron 3 Ultra via OpenRouter (plus four calls that failed: Kimi K2.6, GLM 5.3, DeepSeek V4.1 Flash, Mistral Magistral) with a prompt containing no patient data. Models were used for analysis code, literature triage, critique and drafting; every reference was verified programmatically against PubMed (68/68). All controlled data were deleted from the workstation on 17 October 2026 (see Section 2.1).

## Summary

Mosaic variegated aneuploidy type 1 (MVA1) is caused by biallelic variants in *BUB1B*, which encodes the spindle
assembly checkpoint pseudokinase BUBR1 [1, 2]. We argue that the right unit for drug
repurposing in MVA1 is not the gene but the **allele class**: truncating alleles make no protein, missense alleles in or
near the kinase domain make protein that is turned over fast, and other missense alleles keep normal abundance but
impair checkpoint function [3]. A dosage-raising drug can only help if one allele still makes protein.
We built a citation-locked framework that classifies alleles, scores eleven approved medicines on six referenced
evidence dimensions, ranks them per genotype and tests robustness over 5,000 perturbed weightings - and in this
submission we **instantiated it for the hackathon proband**.

**The proband's genotype.** Local analysis of the provided VCF found exactly two rare protein-altering *BUB1B*
variants, both heterozygous: **c.2210T>G p.(Leu737Ter)**, a ClinVar pathogenic nonsense allele predicted to trigger
nonsense-mediated decay, and **c.3006T>G p.(Asn1002Lys)**, a novel missense change (absent from gnomAD) inside the
pseudokinase domain, ten residues from the HSP90-dependent MVA1 allele L1012P. Phase cannot be resolved from the data
(10.9 kb apart); trans configuration is assumed. No other rare damaging variant was found in 15 further checkpoint/MVA
genes. In the AlphaFold model, Asn1002 is a well-packed (pLDDT 91, 16% exposed) residue whose side chain makes three
hydrogen bonds to backbone atoms of Trp978 and Val998; ESM-2 and ThermoMPNN place N1002K among the published unstable
MVA1 alleles (ESM-2 log-likelihood ratio -5.8, below every ClinVar-benign BUB1B missense; ThermoMPNN +1.5 kcal/mol, as
destabilising as L844F), AlphaMissense calls it likely pathogenic (0.92), while ESM-1v calls it tolerated - an honest
disagreement we report rather than hide. The proband therefore carries the configuration the framework was built for:
**one null allele plus one probably unstable, protein-making kinase-domain allele**.

**Patient-specific ranking.** Treating P(missense allele is unstable) = 0.8 (rule of succession on the 7/8 concordance
of the kinase-window rule) and adding phenotype-driven risk modifiers (rhabdomyosarcoma therapy, nephrocalcinosis), the
leading hypotheses for this child are: **H1** NAD+ repletion with niacin to support SIRT2-dependent stabilisation of the
N1002K protein (rank 1 in 72% and top three in 79% of 5,000 perturbed rankings that also resample the missense class);
**H2** heat-shock-response amplification with arimoclomol to help the N1002K protein fold (rank 2 if instability is
confirmed); **H3** low-dose mTORC1 dampening (everolimus) or IL-1 blockade (anakinra) against progeroid, SASP-driven
decline; **H4** intermittent senolytic clearance (dasatinib). The whole ranking hinges on one measurable quantity - the
half-life of BUBR1-N1002K: niacin leads whenever P(unstable) > 0.5, and if the protein proves stable, allele-independent
strategies take over. We give a decisive, cheap experiment that settles this, a combination rationale (niacin plus
arimoclomol act on two different degradation routes of the same allele), a patient-specific **do-no-harm layer**
(spindle poisons in standard rhabdomyosarcoma regimens, HSP90 inhibitors that would strip the N1002K protein of its
chaperone, nephrotoxic and CYP3A-interacting agents) and an ideation-and-critique round with four other AI models
whose adopted suggestions were all re-verified against PubMed. All references are verified programmatically at build
time; the build fails on any unverifiable citation.

## 1. Variant mechanism: what goes wrong, and where a medicine could act

**Loss of function by dosage.** MVA1 is recessive and biallelic; the shared consequence of the two alleles is too
little functional BUBR1 [1, 3]. In biallelic patient cell lines, a missense allele usually
pairs with a truncating allele, and the cells show an impaired mitotic checkpoint, chromosome alignment defects and
low overall BUBR1; re-expressing BUBR1 restores checkpoint activity, which establishes BUBR1 dysfunction as causal
[3]. We characterise the disease as **loss of function with a dosage component**; in the sources we
reviewed, no evidence supports a gain-of-function mechanism for MVA1 alleles.

**Two allele classes, two kinds of missense.** Truncating alleles produce no transcript; missense alleles in or near
the kinase domain produce protein with high turnover; other missense alleles keep normal abundance but impair
checkpoint activity or microtubule attachment [3]. Selected missense mutants (I909T, L1012P) depend
on HSP90 to escape proteasomal degradation, i.e. they behave as folding-defective proteins [3].
Regulatory alleles that lower expression of wild-type protein have also been described, including an intergenic
variant upstream of *BUB1B* in Japanese families [4, 5].

**Pathway and downstream consequences** (Figure 1; evidence level in brackets: 3 = shown in BubR1-insufficient mice or
MVA patient cells, 2 = shown in cells or one model genotype, 1 = inference):

1. Low BUBR1 -> weak spindle checkpoint and attachment errors [3] [3, 6].
2. -> chromosome missegregation and mosaic aneuploidy [3] [1, 7].
3. -> p53-p21 and p16Ink4a senescence with a senescence-associated secretory phenotype (SASP) [3]
   [8, 9, 10].
4. -> progeroid tissue decline: cataract, lipodystrophy, sarcopenia [3] [7, 11].
5. -> embryonal cancers (rhabdomyosarcoma, Wilms tumour, leukaemia) [3] [1, 12].
6. Allele-specific branch: heterozygous mice carrying the truncating X753 allele show skeletal-muscle mTORC1
   hyperactivity (raised p70 S6 kinase and 4EBP1 phosphorylation) that tracks with sarcopenia, while carriers of a
   missense allele with a similar BUBR1 reduction do not [2] [13].
7. Counter-regulation: NAD+-dependent SIRT2 deacetylation maintains BUBR1 levels; NAD+ decline lowers BUBR1, and the
   NAD+ precursor NMN raised BUBR1 abundance in vivo [3] [14, 15].
8. Aneuploid, arrested cells can be cleared by NK cells via an NF-kB-driven program [2] [10, 16].

The genetics are broader than *BUB1B*: biallelic *CEP57*, *TRIP13*, *BUB1* and *MAD1L1* variants cause related
aneuploidy syndromes [17, 18, 19, 20, 21]. Clinically, MVA
is heterogeneous, with distinct subtypes and cancer predisposition [22, 23].
Increased BubR1 expression protects mice against aneuploidy and cancer and extends healthy lifespan [6],
which makes **raising residual BUBR1** the most causal lever and **limiting downstream senescence/SASP** the most
allele-independent one.

![Figure 1. Mechanism map: the MVA1 disease chain and the intervention points of strategies S1-S3 (S4 acts on aneuploid clones). Solid arrows: level-3 evidence; dashed: one mouse genotype. The do-no-harm layer is shown at the bottom.](https://github.com/ttimesai-star/mva-allele-aware-repurposing/raw/main/figures/mechanism_map.png)

## 2. Patient-specific analysis

### 2.1 Data handling

The controlled files used were the single-sample VCF (GRCh38, Sentieon/GATK) with its index, the clinical phenotype
reference document (HPO table) and the dataset README; the FASTQ files were not downloaded. All files were held on one
separate workstation, never on the authoring laptop, in a git repository or in any chat, and were processed by local
scripts (`src/vcf_panel_extract.py`) that stream the VCF once and keep only records in exons (+/-20 bp) and 2 kb
promoters of 16 spindle-checkpoint/MVA genes. Those 102 panel records - not the genome - were annotated with the
Ensembl Variant Effect Predictor REST service [24], which the organisers' guidance classes as a processing
tool. What left the workstation was a derived summary: four rare protein-altering panel calls (the two *BUB1B* variants
below and two benign-predicted missense variants in other genes, not named here) with zygosity and public annotations,
per-gene variant counts, and the HPO terms. Phase was never resolvable and no secondary-findings search was run. The
VCF, index, phenotype document, intermediate tables and caches were deleted after the analysis (17 October 2026, 21:03 UTC+3; a disk-wide search afterwards found no remaining copy); the
deletion will be confirmed to Sage Bionetworks by e-mail as the rules require. Sequence and structure modelling
(Section 2.3) used only public inputs (UniProt sequence, AlphaFold model) plus the protein change N1002K.

### 2.2 Genotype

| | Allele 1 | Allele 2 |
|---|---|---|
| Change (NM_001211.6) | c.2210T>G p.(Leu737Ter) | c.3006T>G p.(Asn1002Lys) |
| Zygosity, quality | heterozygous, PASS, GQ 99 | heterozygous, PASS, GQ 99 |
| Location | exon 17 of 23, before the pseudokinase domain | exon 23 (last), pseudokinase domain C-lobe |
| Population frequency | gnomAD max AF 7.9e-5 | absent from gnomAD |
| ClinVar | Pathogenic/Likely pathogenic, multiple submitters, MVA1 | not listed; same protein change via c.3006T>A is a VUS |
| Predicted consequence | UGA stop followed by A; nonsense-mediated decay expected | AlphaMissense 0.92 (likely pathogenic); REVEL 0.47; CADD 24.5 |
| Framework class | TRUNC_NONSENSE (no protein) | MS_UNSTABLE (probable fast turnover) |

Table 2. The proband's *BUB1B* alleles (`config/patient_genotype.yaml`; annotations from VEP with AlphaMissense,
REVEL, CADD and SpliceAI [25, 26, 27, 28]).

The nonsense allele sits far upstream of the last exon junction, so its transcript should be degraded by
nonsense-mediated decay [29]; even an escaping product would lack the entire pseudokinase domain. The two
variants are 10.9 kb apart, beyond short-read phasing, there are no parental data, and the organisers declined to
disclose phase during the challenge; compound heterozygosity is assumed from the recessive MVA1 phenotype. Under ACMG
criteria [30] the missense allele has PM2 (absent from population data), PP3 (computational evidence; see
below) and, conditional on trans phase, PM3; it remains formally a variant of uncertain significance - exactly the
"invisible" allele class described in Section 4.1. The panel screen found no other rare damaging variant in 15 further
spindle-checkpoint/MVA genes (*CEP57*, *TRIP13*, *CENATAC*, *BUB1*, *BUB3*, *MAD1L1*, *MAD2L1*, *TTK*, *CDC20*, *KNL1*,
*ZWINT*, *CENPE*, *SGO1*, *KIF2C*, *ANAPC1*).

![Figure 2. The proband's alleles (above the bar) against published MVA1 missense alleles (below; dark = low protein abundance in patient or transfected cells). N1002K lies ten residues from L1012P, which is HSP90-dependent and degraded by the proteasome.](https://github.com/ttimesai-star/mva-allele-aware-repurposing/raw/main/figures/patient_alleles.png)

### 2.3 Is the N1002K protein unstable? Structure and sequence evidence

This is the question the ranking depends on, because a dosage- or folding-directed drug needs a protein-making allele
whose product is degraded too fast.

* **Position.** The pseudokinase domain of human BUBR1 is needed for conformational stability rather than catalysis
  [2]. Four of the five published low-abundance MVA1 alleles lie in it, including L1012P ten residues
  away [3]; mice carrying the mouse equivalent of L1012P (BubR1 L1002P in mouse numbering) together
  with a truncating allele - the proband's configuration - died prematurely [13].
* **Structure.** In the AlphaFold model of BUBR1 [31, 32] Asn1002 has high confidence (pLDDT 91),
  is largely buried (relative solvent accessibility 0.16, 18 C-alpha neighbours within 8 A) and its side-chain amide
  makes three hydrogen bonds with backbone atoms of Trp978 and Val998. Lysine cannot reproduce the acceptor bond from
  OD1 to the Trp978 amide and brings a long charged side chain into a packed site - a structural rationale for
  destabilisation, not proof of it.
* **Stability prediction.** ThermoMPNN [33] on the kinase-domain crop (residues 720-1050) predicts
  +1.5 kcal/mol for N1002K (positive = destabilising in this output), more than for 73% of all possible substitutions in
  the domain, equal to L844F (+1.5), below I909T (+1.9) and L1012P (+3.1), and far above the four ClinVar-benign
  kinase-domain variants (+0.08 to +0.22).
* **Protein language models.** ESM-2 [34] gives N1002K a masked-marginal log-likelihood ratio of -5.8, lower than
  every one of 28 ClinVar benign/likely benign *BUB1B* missense variants (lowest -3.4) and within the range of the
  published unstable alleles (-2.8 to -6.9). ESM-1v (one ensemble member; model family of [35]) disagrees:
  +0.59, inside the benign range, although it does score I909T (-5.4) as damaging. AlphaMissense calls N1002K likely
  pathogenic (0.92).

Three of four predictors and the structural context support instability; one predictor does not. None of these tools
was calibrated on BUBR1 turnover, so we do **not** let them move the prior: P(unstable) stays at 0.8 from the
kinase-window rule (7 of 8 published alleles concordant; (7+1)/(8+2) by Laplace's rule of succession), and the class
question is carried through the ranking as uncertainty (Section 2.4) and handed to an experiment (Section 2.5).

![Figure 3. Sequence and structure models place N1002K (red line) with the published low-abundance MVA1 alleles in ESM-2 and ThermoMPNN, but with benign variants in ESM-1v (outputs/variant_effect_scores.csv).](https://github.com/ttimesai-star/mva-allele-aware-repurposing/raw/main/figures/variant_effect.png)

### 2.4 Ranking for this child

The patient-specific score mixes the two possible classes of the missense allele
(score = 0.8 x score(nonsense + unstable missense) + 0.2 x score(nonsense + stable missense)) and adds phenotype-driven
risk modifiers from the HPO profile: +1.0 risk for gentamicin because of nephrocalcinosis (HP:0000121; aminoglycoside
nephrotoxicity affects 20-33% of exposed children [36]); +0.5 for dasatinib (myelosuppression on top of
rhabdomyosarcoma chemotherapy, HP:0002859, which already required dose reduction in an MVA1 child
[37]); +0.5 for everolimus, sirolimus, ruxolitinib and anakinra (immunosuppression after cancer therapy
and possible loss of immune clearance of aneuploid cells [16]).

| Rank | Candidate (strategy) | If N1002K unstable | If N1002K stable | Expected (P = 0.8) | P(rank 1) | P(top 3) |
|---|---|---|---|---|---|---|
| 1 | niacin (S1, BUBR1 dosage) | 0.67 | 0.37 | **0.61** | **72%** | **79%** |
| 2 | everolimus (S2, mTORC1/SASP) | 0.52 | 0.52 | 0.52 | 10% | 51% |
| 3 | dasatinib (S3, senolytic) | 0.48 | 0.48 | 0.48 | 5% | 36% |
| 4 | metformin (S4, aneuploid stress) | 0.47 | 0.47 | 0.47 | 6% | 24% |
| 5 | sirolimus (S2, mTORC1/SASP) | 0.47 | 0.47 | 0.47 | 2% | 21% |
| 6 | arimoclomol (S1b, folding) | 0.53 | 0.17 | 0.46 | 2% | 46% |
| 7 | hydroxychloroquine (S4, aneuploid stress) | 0.45 | 0.45 | 0.45 | 1% | 13% |
| 8 | anakinra (S2b, IL-1/SASP) | 0.43 | 0.43 | 0.43 | 1% | 9% |
| 9 | sodium phenylbutyrate (S1b, folding) | 0.49 | 0.16 | 0.42 | 1% | 19% |
| 10 | gentamicin (S1c, read-through) | 0.38 | 0.38 | 0.38 | 0% | 2% |
| 11 | ruxolitinib (S2, mTORC1/SASP) | 0.33 | 0.33 | 0.33 | 0% | 0% |

Table 3. Patient-specific ranking (`outputs/patient_ranking.csv`) and stability over 5,000 draws that perturb weights
and expert scores and, in every draw, resample the class of the missense allele (`outputs/patient_rank_stability.csv`).

Three results matter. First, niacin is the most robust single hypothesis for this child, but only because the
N1002K allele probably makes protein; it ranks first whenever P(unstable) exceeds 0.5 (Figure 4). Second, arimoclomol
has the widest spread: second if instability is confirmed, near the bottom if not - it is the candidate whose value an
experiment would change most. Third, removing the phenotype modifiers changes little at the top
(`outputs/patient_ranking_no_context.csv`); they mainly push gentamicin, already low, further down and separate the
immunosuppressive senomorphics from metformin.

![Figure 4. Expected score of the leading candidates as a function of P(N1002K protein is unstable). Niacin overtakes everolimus at P = 0.5; arimoclomol overtakes everolimus only near certainty.](https://github.com/ttimesai-star/mva-allele-aware-repurposing/raw/main/figures/patient_class_sweep.png)

### 2.5 The decisive experiment and a combination rationale

**One assay decides the ranking:** the half-life of BUBR1-N1002K. Cheapest version, no patient material needed:
express wild-type, N1002K and L1012P (positive control) BUBR1 in a human cell line, measure turnover by cycloheximide
chase, and test HSP90 dependence with an HSP90 inhibitor, as was done for L1012P [3]. If N1002K behaves
like L1012P, test niacin (or another NAD+ precursor), arimoclomol and their combination for steady-state BUBR1 and
checkpoint competence (nocodazole arrest, mitotic-error rate). If the treating team can provide patient fibroblasts or
lymphoblastoid cells, the same readouts on endogenous protein are definitive.

**Why a combination is worth testing.** An external reviewer model objected that SIRT2-dependent deacetylation, the
target of the NAD+ route, controls degradation of correctly folded BUBR1 [14, 15] and cannot rescue
a folding defect. We agree this is the main weakness of H1 for a misfolded allele. Niacin and arimoclomol act on two
different degradation routes of the same protein - acetylation-dependent turnover and folding quality control
(arimoclomol co-induces heat-shock proteins [38]) - so the combination, tested in the assay above, is the
patient-specific hypothesis with the best mechanistic coverage. Combinations of approved medicines are in scope for
Track 2.

**Nonsense allele.** The UGA stop of p.(Leu737Ter) is followed by A; aminoglycoside read-through depends on the stop
codon and its flanking context [39] and on transcripts escaping decay. Gentamicin therefore stays a
cell-assay probe only, now further penalised for this child's nephrocalcinosis.

### 2.6 Patient-specific do-no-harm notes

These are observations for the treating oncology team to weigh, not recommendations to withhold standard care.

1. **Spindle poisons in rhabdomyosarcoma regimens.** Standard regimens contain vincristine, a tubulin binder;
   BubR1-insufficient cells slip out of mitotic arrest [40], and an MVA1 child with rhabdomyosarcoma was managed
   with reduced-intensity chemotherapy [37].
2. **HSP90 inhibitors.** If N1002K behaves like L1012P, HSP90 inhibition would strip it of its chaperone and lower the
   only protein-making allele [3].
3. **ATP-site kinase inhibitors that bind BUBR1.** ATP-competitive inhibitors can deprive Hsp90-dependent kinases of
   Cdc37 and drive their degradation [41]. ChEMBL lists no approved drug with meaningful BUBR1 affinity
   (24 activities; best approved ceritinib, IC50 13 uM), so this is a watch item, not a current conflict.
4. **Renal.** Nephrocalcinosis plus aminoglycosides: avoid gentamicin read-through outside cell assays
   [36]; renal ultrasound surveillance is already recommended in MVA [12].
5. **Drug interactions.** Vincristine is CYP3A-metabolised and strong CYP3A inhibitors cause severe vincristine
   neurotoxicity [42]; mTOR inhibitors and dasatinib need an interaction review against the regimen.

### 2.7 Multi-model ideation and critique

To reduce single-model bias we sent a structured review request - containing only a generic genotype class
(nonsense + kinase-domain missense) and our candidate list, no patient data - to six other models through free tiers;
every call, including failures, is logged verbatim in `logs/multimodel_critique/`. Two returned reviews (Google Gemini
3.7 Flash; NVIDIA Nemotron 3 Ultra via OpenRouter); four failed (Moonshot Kimi K2.6 and GLM 5.3 returned empty
content, DeepSeek V4.1 Flash timed out, Mistral Magistral was rate-limited). Both reviewers independently made the same
main objection: SIRT2-dependent deacetylation, the target of the NAD+ route, controls degradation of correctly folded
BUBR1 [14] and need not rescue a folding defect - which we accept and answer with the combination rationale in
Section 2.5. Adopted after our own PubMed verification: IL-1 receptor blockade with anakinra as a lower-toxicity
senomorphic (new candidate H3b) [43, 44, 45]; the vincristine interaction flag
[42]; and a new counter-argument to sodium phenylbutyrate - inhibiting HDAC6 hyperacetylates HSP90 and
releases its clients [46], which could work against an HSP90-dependent N1002K protein (whether phenylbutyrate
does this at clinical exposure is unknown). Not adopted: suggestions without verifiable evidence for a BUBR1-relevant
mechanism or with unacceptable risk in a cancer-predisposed child (a proteasome inhibitor, ursodeoxycholic acid,
an eIF2-alpha phosphatase inhibitor, ACE inhibitors and statins as senomorphics); ataluren, which is not
FDA-approved and duplicates the read-through probe; and the reviewers' citations, which we did not take on trust - one
PMID given for the ataluren discovery paper pointed to an unrelated article (the correct record is PMID 17450125).

### 2.8 A discovery lead the genotype suggests (not repurposing)

Because BUBR1's ATP pocket serves stability rather than catalysis [2], a pocket ligand would not cost
function and might thermodynamically stabilise N1002K; but the same kind of ligand can also remove kinases from
Hsp90-Cdc37 and promote their degradation [41]. No approved drug binds BUBR1 well enough to test this by
repurposing; we record it as a two-edged question for structure-based design (`config/evidence.yaml`,
`discovery_directions`).

## 3. Approach

The pipeline (`run_all.py`) is fully scripted and reruns offline from a committed HTTP cache.

1. **Reference verification** (`verify_refs.py`). Each of 68 references is resolved against PubMed E-utilities; title
   word overlap and DOI must match. Any citation key used in the knowledge base but not defined, or any reference that
   fails verification, stops the build. Result: 68/68 verified (`outputs/reference_verification.csv`).
2. **Allele-class landscape** (`clinvar_landscape.py`). All single-gene Pathogenic/Likely pathogenic ClinVar records for
   seven MVA genes are classified as truncating (frameshift/splice), nonsense, kinase-window missense (unstable),
   other missense (functional) or regulatory, using the UniProt O60566 kinase domain (766-1050) extended N-terminally to
   residue 720 so that R727C, which is unstable in patient cells, counts as "near" [3].
3. **Rule check against measured protein abundance.** The kinase-window rule is compared with the eight MVA1 missense
   alleles whose protein abundance was measured [3].
4. **Drug annotation** (`drug_annotations.py`). Each candidate is checked live in ChEMBL (max phase 4, first approval,
   mechanism) and openFDA (label section 8.4 "Pediatric Use"); Open Targets supplies MVA-associated genes for scaling.
5. **Do-no-harm scan.** Every approved (ChEMBL max phase 4) medicine whose mechanism matches a checkpoint-dependent
   pattern (tubulin, Aurora kinase, Polo-like kinase, KIF11/Eg5, HSP90) is listed automatically.
6. **Scoring** (`score.py`). Each candidate gets expert scores (0-3) for causal proximity to the BUBR1 defect (*mech*),
   evidence in BubR1/MVA models (*model*), human target-engagement or efficacy evidence (*human*), paediatric
   experience and safety (*ped*), speed and cost of a decisive test (*test*) and a specific risk for an MVA child
   (*risk*, subtracted). Every score carries references in `config/evidence.yaml`.
   `base = sum(w_i * s_i / 3) - w_risk * risk / 3`, with weights mech 0.25, model 0.20, human 0.15, ped 0.20,
   test 0.10, risk 0.10. Genotype fit is the best fit of the strategy to either allele (one responsive allele is enough
   for a dosage effect); a mismatched strategy keeps 25% of its score: `score = base * (0.25 + 0.75 * fit)`.
7. **Robustness.** 5,000 draws with Dirichlet-perturbed weights (concentration 20) and uniform noise of +/-0.5 on every
   expert score; we report P(rank 1), P(top 3) and median rank.
8. **Genotype instantiation** (`classify_variants.py`). Given two alleles, the script prints their classes and the
   matching ranked list.
9. **Patient layer** (Submission 2). `vcf_panel_extract.py` (run where the controlled VCF lives) finds rare
   protein-altering variants in 16 checkpoint/MVA genes; `variant_effect_models.py` scores missense variants with
   ESM-1v, ESM-2 and ThermoMPNN on the public AlphaFold model; `score.py` computes the class-mixture ranking with
   phenotype risk modifiers, a 5,000-draw stability analysis that also resamples the missense class, and a sweep
   over P(unstable); `patient_figures.py` draws Figures 2-4.

## 4. Results

### 4.1 The dosage-responsive allele class is invisible in ClinVar

Of 78 single-gene Pathogenic/Likely pathogenic *BUB1B* records, 76 are truncating (40 frameshift/splice, 36 nonsense),
one is a missense outside the kinase window and one is unclassified (`outputs/clinvar_allele_class_summary.csv`). In
contrast, of 421 *BUB1B* missense records linked to MVA or premature chromatid separation, 392 are of uncertain
significance and 22 have conflicting classifications; one is Likely pathogenic
(`outputs/bub1b_missense_clinvar_classes.csv`). The missense half of a typical MVA1 genotype - exactly the half that
determines whether a dosage-raising drug can work - is therefore usually a variant of uncertain significance in public
databases. The kinase-window rule reproduces measured protein abundance for 7 of 8 published MVA1 missense alleles;
the single discordant allele (Q921H, normal abundance) occurs in cis with the unstable L844F allele
[3] (`outputs/kinase_window_rule_check.csv`, Figure 5).

![Figure 5. BUBR1 allele map. Above: ClinVar Pathogenic/Likely pathogenic variants with a protein position. Below: published MVA1 missense alleles coloured by measured protein abundance; the shaded window marks the "in or near kinase domain" rule.](https://github.com/ttimesai-star/mva-allele-aware-repurposing/raw/main/figures/bubr1_alleles.png)

### 4.2 The ranking depends on the genotype

| Candidate (strategy) | trunc + unstable missense | nonsense + unstable missense | trunc + stable missense | trunc + low-expression | two unstable missense |
|---|---|---|---|---|---|
| niacin (S1) | **0.67** | **0.67** | 0.37 | **0.57** | **0.67** |
| arimoclomol (S1b) | 0.53 | 0.53 | 0.17 | 0.13 | 0.53 |
| everolimus (S2) | 0.53 | 0.53 | **0.53** | 0.53 | 0.37 |
| dasatinib (S3) | 0.50 | 0.50 | 0.50 | 0.50 | 0.50 |
| sodium phenylbutyrate (S1b) | 0.49 | 0.49 | 0.16 | 0.12 | 0.49 |
| sirolimus (S2) | 0.48 | 0.48 | 0.48 | 0.48 | 0.34 |
| metformin (S4) | 0.47 | 0.47 | 0.47 | 0.47 | 0.47 |
| hydroxychloroquine (S4) | 0.45 | 0.45 | 0.45 | 0.45 | 0.45 |
| ruxolitinib (S2) | 0.35 | 0.35 | 0.35 | 0.35 | 0.24 |
| anakinra (S2b) | 0.45 | 0.45 | 0.45 | 0.45 | 0.45 |
| gentamicin (S1c) | 0.10 | 0.42 | 0.10 | 0.10 | 0.10 |

Table 4. Allele-aware scores for generic genotype configurations (`outputs/ranking_by_genotype.csv`; Figure 6). The first column is the most common MVA1
configuration [3] and our primary scenario. Scores rank hypotheses for testing; they are not
predictions of efficacy.

![Figure 6. Allele-aware score by BUB1B genotype configuration.](https://github.com/ttimesai-star/mva-allele-aware-repurposing/raw/main/figures/ranking_heatmap.png)

### 4.3 Robustness

Under 5,000 perturbations of weights and expert scores, niacin ranked first in 87.8% and in the top three in 97.4% of
draws for the primary generic genotype. Arimoclomol (52.4% top three) and everolimus (50.5%) compete for second place,
followed by dasatinib (33.4%); ruxolitinib and gentamicin never entered the top three
(`outputs/rank_stability_G1_TRUNC_x_MSUNSTABLE.csv`, Figure 7). The top of the list is robust to reasonable
disagreement about weights; places two to four are not, and we present them as a group.

![Figure 7. Rank stability for the primary genotype across 5,000 perturbed rankings.](https://github.com/ttimesai-star/mva-allele-aware-repurposing/raw/main/figures/rank_stability.png)

## 5. Repurposing hypotheses

Evidence tier for each hypothesis: **A** - the intervention itself was tested in a BubR1/MVA model; **B** - the
mechanism was validated (genetically or with a related agent) in a BubR1/MVA model and the drug engages that
mechanism in humans; **C** - mechanistic inference from other systems; **D** - exploratory. None of the hypotheses
below has tier A evidence.

### H1. NAD+ repletion with niacin to stabilise residual BUBR1 (tier B; strategy S1)

*Rationale.* BUBR1 declines with NAD+; SIRT2 keeps BUBR1 deacetylated and abundant; the NAD+ precursor NMN raised
BUBR1 abundance in vivo, and SIRT2 overexpression extended median lifespan of BubR1 hypomorphic mice [14];
SIRT2 deacetylation also governs BUBR1 degradation during mitosis [15]. More BUBR1 protects against
aneuploidy and cancer [6]. Niacin is a licensed NAD+ precursor that raised blood NAD+ up to 8-fold and
restored muscle NAD+ to control levels in adult patients with mitochondrial myopathy [47] and was used in a
child with NAD+ deficiency [48]. We prefer niacin to nicotinamide because nicotinamide is itself a sirtuin
inhibitor [49].

*Genotype dependence.* Only meaningful if at least one allele makes protein (unstable missense, regulatory or
functional missense); predicted to be inert for two null alleles.

*Counter-evidence and gaps.* Niacin itself has never been tested on BUBR1; SIRT2 overexpression did not extend lifespan
in wild-type mice [50]; paediatric niacin experience (hypercholesterolaemia) showed frequent reversible adverse
effects, including aminotransferase elevations [51], and the US label states that paediatric safety is not
established. Raising NAD+ in a cancer-predisposed child has theoretical tumour-metabolism implications.

*Decisive test.* BUBR1 western blot, nocodazole checkpoint assay and mitotic-error rate in MVA1 patient fibroblasts or
lymphoblastoid cells with and without an NAD+ precursor, stratified by allele class.

### H2. Heat-shock-response amplification with arimoclomol to rescue a folding-defective missense BUBR1 (tier C; S1b)

*Rationale.* Kinase-domain missense mutants turn over about twice as fast as wild type in patient cells and depend on
HSP90 to avoid proteasomal clearance [3]. Arimoclomol amplifies the heat-shock response to counter
protein misfolding; in a randomised placebo-controlled trial in patients aged 2-18 years with Niemann-Pick disease
type C it reduced annual disease progression by 65% [52], and heat-shock protein amplification improved
cerebellar myelination in an NPC mouse model [53]. Sodium phenylbutyrate, a licensed chemical chaperone, is a
weaker alternative: its rescue evidence is for a secretory protein (LDLR G544V) [54], whereas BUBR1 is
cytosolic, and it also inhibits histone deacetylases, which could alter BUBR1 acetylation [55].

*Genotype dependence.* Only for unstable missense alleles; near-zero fit otherwise (Table 4).

*Counter-evidence and gaps.* No data on BUBR1; heat-shock proteins also support tumour cells (theoretical);
treatment-related urticaria and angioedema were reported [52].

*Decisive test.* Cycloheximide chase of BUBR1 in patient cells carrying a kinase-domain missense allele, with and
without arimoclomol; readout BUBR1 half-life and checkpoint competence.

### H3. Low-dose intermittent mTORC1 dampening (everolimus; sirolimus as alternative) against progeroid decline (tier B-/C; S2)

*Rationale.* Truncating-allele (X753) carrier mice show skeletal-muscle mTORC1 hyperactivity that tracks with
sarcopenia [13]; mTOR drives the pro-tumorigenic SASP via IL1A translation [43]. Low-dose
everolimus regimens improved, rather than suppressed, influenza-vaccine responses and reduced infections in older
adults [56, 57]. Everolimus carries a paediatric label from age 1 (TSC-associated SEGA; openFDA).

*Scope.* Proposed for tissue ageing (muscle, fat), not as an anticancer drug: temsirolimus had activity in relapsed
rhabdomyosarcoma [58] but did not improve event-free survival when added to frontline therapy for
intermediate-risk rhabdomyosarcoma [59].

*Counter-evidence and gaps.* Rapamycin analogues have not been tested in BubR1 mice in the sources reviewed;
immunosuppression could blunt NK-cell clearance of aneuploid cells [16].

*Decisive test.* p-S6K/p-4EBP1 and SASP panel in MVA fibroblasts; muscle function in BubR1 +/X753 mice.

### H3b. IL-1 receptor blockade (anakinra) as a lower-toxicity senomorphic (tier C; S2b; added in Submission 2)

*Rationale.* mTOR drives the SASP by promoting IL1A translation [43], so blocking IL-1 signalling targets the
same output as H3 one step further downstream, without mTOR inhibition. Anakinra is used from the neonatal period in
NOMID with sustained benefit [44, 45]. Suggested in the multi-model critique round and kept
after verification. *Counter-evidence and gaps.* No BubR1/MVA data; the SASP also recruits immune clearance of
aneuploid senescent cells [16]. *Decisive test.* IL-1a/IL-6/IL-8 secretion and SA-beta-gal in MVA fibroblast
cultures with and without IL-1 receptor blockade.

### H4. Intermittent senolytic clearance (dasatinib, with quercetin) of senescent aneuploid cells (tier B; S3)

*Rationale.* In BubR1 progeroid mice, genetic removal of p16Ink4a-positive cells delayed cataract, lipodystrophy and
sarcopenia, and late-life clearance attenuated established disease [11]; p16 and p21 shape these phenotypes
[8, 9]. Dasatinib plus quercetin is the senolytic pair identified by transcriptome-guided
screening [60] and reduced senescent-cell burden in a first human study [61]. Dasatinib has
paediatric experience in chronic myeloid leukaemia [62]. This strategy is allele-independent.

*Counter-evidence and gaps.* Senolytic drugs were not tested in BubR1 mice in our sources; dasatinib causes
myelosuppression and effusions; quercetin is a supplement, not a licensed medicine (fisetin is another non-licensed
senolytic [63]). Arrested aneuploid cells may already be cleared by the immune system
[10, 16], which lowers the expected added benefit.

*Decisive test.* SA-beta-gal and p16 counts in MVA fibroblast cultures after pulsed exposure; cataract and fat
endpoints in BubR1 H/H mice.

### H5. Exploratory hypotheses (tier C-D)

* **Aneuploid-cell stress (S4).** Chloroquine (autophagy inhibition) and AICAR (energy stress) were selectively
  anti-proliferative in trisomic mouse cells and aneuploid cancer lines [64]. Hydroxychloroquine and metformin
  (an AMPK activator; inference only) are licensed proxies. No MVA data.
* **Read-through of a nonsense allele (S1c).** Intravenous gentamicin produced read-through in children with nonsense
  variants in another disease [65]. Relevant only if one allele is a premature stop codon (36 of 78
  ClinVar P/LP *BUB1B* records); nephro- and ototoxicity make this a cell-assay probe, not a long-term option.
* **Ruxolitinib.** JAK inhibition suppressed the SASP and frailty in old mice [66]; weakest link to BUBR1 and a
  high infection risk; ranked last but one.

## 6. Do-no-harm layer: risks for the clinician to weigh

These flags are generated from mechanism and literature. They are research observations, not prescribing advice,
and are **not** a recommendation to withhold standard-of-care oncology treatment, which remains the decision of the
treating oncology team.

1. **Spindle-poison antimitotics.** BubR1-insufficient cells are defective in checkpoint activation and slip out of
   mitotic arrest [40]. The automated scan lists 23 approved tubulin-targeting agents in ChEMBL
   (`outputs/liability_scan.csv`), including vinca alkaloids and taxanes used in paediatric oncology, as well as
   benzimidazole anthelmintics (albendazole, mebendazole, thiabendazole) and griseofulvin, which bind tubulin.
   Whether the latter matter clinically in MVA is unknown; we flag them because they are commonly used in children.
2. **Severe treatment toxicity.** Case reports describe MVA1 children with Wilms tumour and with rhabdomyosarcoma, the
   latter managed with reduced-intensity chemotherapy [37, 67]; genomic instability disorders
   can show exquisite sensitivity to chemotherapy and ionising radiation [12].
3. **Radiation in surveillance.** The AACR workshop update supports renal ultrasound surveillance every 3 months from
   birth until age 7 for all MVA conditions and avoidance of radiation exposure [12].
4. **HSP90 inhibitors.** Kinase-domain missense BUBR1 needs HSP90 to escape degradation [3]; HSP90
   inhibition would be expected to lower it further. Pimitespib is approved in Japan [68]. An HSP90 inhibitor
   (17-AAG) was also among the aneuploidy-selective compounds [64], so this is a genuine trade-off for the
   clinician, not a settled contraindication.
5. **Mitotic kinase and Eg5 inhibitors** (Aurora, Polo-like kinase, KIF11): no approved agents matched in ChEMBL at the
   time of the scan; flagged for investigational use.
6. **Candidate-specific risks.** Immunosuppression with mTOR/JAK inhibitors (could blunt NK clearance of aneuploid
   cells [16]); myelosuppression with dasatinib; hepatotoxicity and flushing with niacin [51];
   nephro- and ototoxicity with gentamicin [65]; retinal toxicity with long-term hydroxychloroquine
   (label).
7. **Cancer predisposition.** Any agent that supports cell survival (NAD+ repletion, heat-shock proteins) must be
   weighed against tumour risk in a cancer-predisposition syndrome before any clinical use.

## 7. Limitations

* **Genotype limits.** Phase of the two *BUB1B* variants is assumed, not shown (10.9 kb apart; no parental data).
  The VCF contains small variants only; structural or copy-number variants and deep-intronic splice variants were not
  assessed, and only 16 checkpoint/MVA genes were screened.
* **The class of N1002K is predicted, not measured.** Three of four predictors and the structural context favour
  instability, ESM-1v does not, and none is calibrated on BUBR1 turnover; the ranking is therefore reported as a
  mixture with an explicit sweep, and the half-life assay in Section 2.5 should replace the prior.
* **Phenotype modifiers are coarse.** We used two HPO terms only and no treatment history (not provided); the treating
  team's current regimen determines real interactions.
* **Expert scores are judgements.** Every score is referenced and the sensitivity analysis shows which conclusions
  survive disagreement, but the scores were assigned by one participant with AI assistance and need independent review.
* **The kinase-window rule rests on eight measured alleles** from one study [3]; a functional assay
  (cycloheximide chase) should replace the rule for any new allele, including N1002K.
* **Literature bias.** The candidate list comes from mechanism-first literature review, not an unbiased screen; drugs
  acting through mechanisms outside the four strategies would be missed.
* **No efficacy data in MVA for any candidate.** Mouse evidence (BubR1 hypomorphs, X753 and L1002P alleles) may not
  translate; some benefits were sex-specific [14] or not reproduced in wild-type mice [50].
* **Database snapshot.** ClinVar, ChEMBL, openFDA and Open Targets were queried in October 2026; the HTTP cache makes
  the snapshot reproducible but not current.

## 8. Scalability

The framework separates generic machinery (reference verification, ClinVar allele-class landscape, ChEMBL/openFDA
annotation, liability scan, scoring, sensitivity analysis) from a disease-specific YAML knowledge base. It already
profiles six further MVA genes (*TRIP13*, *CEP57*, *BUB1*, *MAD1L1*, *CENATAC*, *BUB3*), and Open Targets lists 11
MVA-associated genes (`outputs/mva_gene_associations.csv`). Adding a disease means writing a new `evidence.yaml`;
adding a patient means two command-line arguments. The pipeline runs in minutes on a laptop without a GPU and reruns
offline.

## 9. Reproducibility

`python run_all.py` reproduces every table and figure (offline from the committed cache after the first run);
`pytest` covers allele classification, citation completeness, genotype fit, deterministic ranking, the recorded
classes of the proband's alleles and the dependence of the patient ranking on the missense class. The patient layer
needs the controlled VCF only for `vcf_panel_extract.py` (rerun on 17 October 2026 from the repository version: same
four calls); `variant_effect_models.py` needs only public inputs and about 45 CPU-minutes. Code is
MIT-licensed; submissions are released under CC BY 4.0 per the hackathon rules. Repository:
https://github.com/ttimesai-star/mva-allele-aware-repurposing (private until the hackathon closes, then public).

## 10. Methods description (Track 2 template)

**Team name.** zbignevich1 (individual participant).

**Approach in detail.** Mechanism-first, allele-aware repurposing: (1) characterise the BUB1B disease chain and the
two functional allele classes from primary literature; (2) define four intervention strategies anchored to nodes of
that chain (raise residual BUBR1; help unstable missense BUBR1 fold; dampen mTORC1/SASP; clear senescent cells) plus
exploratory ones; (3) select approved medicines acting on each node and confirm approval status (ChEMBL max phase 4)
and paediatric label text (openFDA); (4) score each candidate on six referenced evidence dimensions; (5) weight by
allele fit for five genotype configurations; (6) test rank stability by Monte Carlo perturbation; (7) scan all
approved drugs for checkpoint-dependent mechanisms to build the do-no-harm layer; (8) Submission 2: extract the
proband's rare protein-altering variants in 16 checkpoint/MVA genes from the provided VCF, classify both *BUB1B*
alleles, assess the missense allele with ESM-1v, ESM-2, AlphaMissense and ThermoMPNN on the AlphaFold model, and
compute a patient-specific ranking that mixes over the missense class and adds HPO-driven risk modifiers;
(9) multi-model critique round, every adopted suggestion re-verified against PubMed.

**Generative AI (required).** AI use: Anthropic Claude Opus 5.5 via Claude Code on a consumer Claude Max subscription, with the 'Help improve our AI models' (training) setting switched OFF on 17 October 2026, before any hackathon data were opened, and no feedback or ratings given on any output - the configuration the organisers' 28 September 2026 guidance accepts for consumer tiers. Raw controlled files (VCF, phenotype document) were processed only by local scripts on a separate workstation; Claude received only derived summaries (the four rare protein-altering calls in a 16-gene panel with zygosity and public annotations, per-gene counts, and the HPO table). Panel exon records (102) were annotated with Ensembl VEP REST, a processing tool. A critique round used free tiers of Google Gemini 3.7 Flash and NVIDIA Nemotron 3 Ultra via OpenRouter (plus four calls that failed: Kimi K2.6, GLM 5.3, DeepSeek V4.1 Flash, Mistral Magistral) with a prompt containing no patient data. Models were used for analysis code, literature triage, critique and drafting; every reference was verified programmatically against PubMed (68/68). All controlled data were deleted from the workstation on 17 October 2026 (see Section 2.1).

**Automated or manual identification.** Hybrid. Candidate selection and expert scores come from manual, AI-assisted
literature review and curation; database annotation, reference verification, the ClinVar landscape, the liability
scan, scoring and sensitivity analysis are automated.

**Manual literature review and curation.** Targeted PubMed searches on BUBR1 biology, BubR1 mouse models, NAD+/SIRT2,
senescence and senolytics, mTOR and SASP, aneuploidy-selective compounds, proteostasis drugs and the paediatric safety
of each candidate. Each mechanism edge, strategy and score in `config/evidence.yaml` carries citation keys;
counter-evidence is recorded explicitly (`counter_refs`).

**Public data only?** No. Hackathon data used: the VCF (16-gene panel extraction only) and the HPO phenotype terms.
FASTQ files were not used.

**Public data sources.** PubMed/E-utilities (literature and reference verification), ClinVar (variant classes),
UniProt O60566 (domain map), ChEMBL (approval status, mechanisms, liability scan), openFDA drug labels (paediatric use),
Open Targets Platform (MVA-associated genes), MONDO (MONDO_0009759), Ensembl VEP REST with gnomAD, AlphaMissense,
REVEL, CADD and SpliceAI; UCSC Genome Browser API (MANE transcripts); AlphaFold DB (AF-O60566-F1, v6); ESM-1v and
ESM-2 weights (Hugging Face); ThermoMPNN code and weights (GitHub, Kuhlman lab).

**Proprietary data.** None.

**Mechanism characterisation.** Loss of function with a dosage component (Section 1): truncating alleles give no
transcript, kinase-domain missense alleles give unstable protein, other missense alleles give specific checkpoint or
attachment defects [3]; pathway: spindle assembly checkpoint -> chromosome missegregation ->
p53/p21/p16 senescence and SASP -> progeroid decline and embryonal cancer, with an allele-specific mTORC1 branch.

**Time and effort.** About 30 person-hours including AI-assisted work; generic pipeline runs in minutes on a laptop;
patient layer ~45 CPU-minutes (ESM, ThermoMPNN) on a desktop CPU; no paid compute.

**Method abstract.** See *Summary*. Strengths: allele-aware rather than gene-level ranking; every claim tied to a
verified reference; explicit counter-evidence; quantified robustness including uncertainty about the patient's
missense class; patient-specific safety layer; one decisive experiment named; scales by configuration.
Limitations: Section 7.

## Acknowledgement

This work was made possible through the Hackathon, organized by Sage Bionetworks in partnership with the MVA Society,
Hugging Face, and BEACON (The Benchmarking, Evaluation, and Assessment Consortium for Science), with prize sponsorship
from AWS and Anthropic. We are deeply grateful to the child and their family who generously contributed their data and
their story to advance research into this rare disease. We acknowledge their trust in making this Hackathon possible.

## References

1. Hanks S et al. Constitutional aneuploidy and cancer predisposition caused by biallelic mutations in BUB1B. *Nat Genet* 2004. PMID 15475955. doi:10.1038/ng1449.
2. Suijkerbuijk SJ et al. The vertebrate mitotic checkpoint protein BUBR1 is an unusual pseudokinase. *Dev Cell* 2012. PMID 22698286. doi:10.1016/j.devcel.2012.03.009.
3. Suijkerbuijk SJ et al. Molecular causes for BUBR1 dysfunction in the human cancer predisposition syndrome mosaic variegated aneuploidy. *Cancer Res* 2010. PMID 20516114. doi:10.1158/0008-5472.CAN-09-4319.
4. Ochiai H et al. TALEN-mediated single-base-pair editing identification of an intergenic mutation upstream of BUB1B as causative of PCS (MVA) syndrome. *Proc Natl Acad Sci U S A* 2014. PMID 24344301. doi:10.1073/pnas.1317008111.
5. Matsuura S et al. Monoallelic BUB1B mutations and defective mitotic-spindle checkpoint in seven families with premature chromatid separation (PCS) syndrome. *Am J Med Genet A* 2006. PMID 16411201. doi:10.1002/ajmg.a.31069.
6. Baker DJ et al. Increased expression of BubR1 protects against aneuploidy and cancer and extends healthy lifespan. *Nat Cell Biol* 2013. PMID 23242215. doi:10.1038/ncb2643.
7. Baker DJ et al. BubR1 insufficiency causes early onset of aging-associated phenotypes and infertility in mice. *Nat Genet* 2004. PMID 15208629. doi:10.1038/ng1382.
8. Baker DJ et al. Opposing roles for p16Ink4a and p19Arf in senescence and ageing caused by BubR1 insufficiency. *Nat Cell Biol* 2008. PMID 18516091. doi:10.1038/ncb1744.
9. Baker DJ et al. p21 both attenuates and drives senescence and aging in BubR1 progeroid mice. *Cell Rep* 2013. PMID 23602569. doi:10.1016/j.celrep.2013.03.028.
10. Santaguida S et al. Chromosome Mis-segregation Generates Cell-Cycle-Arrested Cells with Complex Karyotypes that Are Eliminated by the Immune System. *Dev Cell* 2017. PMID 28633018. doi:10.1016/j.devcel.2017.05.022.
11. Baker DJ et al. Clearance of p16Ink4a-positive senescent cells delays ageing-associated disorders. *Nature* 2011. PMID 22048312. doi:10.1038/nature10600.
12. Nakano Y et al. Update on Recommendations for Cancer Screening and Surveillance in Children with Genomic Instability Disorders. *Clin Cancer Res* 2024. PMID 39264246. doi:10.1158/1078-0432.CCR-24-1098.
13. Sieben CJ et al. BubR1 allelic effects drive phenotypic heterogeneity in mosaic-variegated aneuploidy progeria syndrome. *J Clin Invest* 2020. PMID 31738183. doi:10.1172/JCI126863.
14. North BJ et al. SIRT2 induces the checkpoint kinase BubR1 to increase lifespan. *EMBO J* 2014. PMID 24825348. doi:10.15252/embj.201386907.
15. Suematsu T et al. Deacetylation of the mitotic checkpoint protein BubR1 at lysine 250 by SIRT2 and subsequent effects on BubR1 degradation during the prometaphase/anaphase transition. *Biochem Biophys Res Commun* 2014. PMID 25285631. doi:10.1016/j.bbrc.2014.09.128.
16. Wang RW et al. Aneuploid senescent cells activate NF-κB to promote their immune clearance by NK cells. *EMBO Rep* 2021. PMID 34105235. doi:10.15252/embr.202052032.
17. Aziz K et al. Mosaic-variegated aneuploidy syndrome mutation or haploinsufficiency in Cep57 impairs tumor suppression. *J Clin Invest* 2018. PMID 30035751. doi:10.1172/JCI120316.
18. Yost S et al. Biallelic TRIP13 mutations predispose to Wilms tumor and chromosome missegregation. *Nat Genet* 2017. PMID 28553959. doi:10.1038/ng.3883.
19. Carvalhal S et al. Biallelic BUB1 mutations cause microcephaly, developmental delay, and variable effects on cohesion and chromosome segregation. *Sci Adv* 2022. PMID 35044816. doi:10.1126/sciadv.abk0114.
20. Villarroya-Beltri C et al. Biallelic germline mutations in MAD1L1 induce a syndrome of aneuploidy with high tumor susceptibility. *Sci Adv* 2022. PMID 36322655. doi:10.1126/sciadv.abq5914.
21. Malumbres M et al. Mosaic variegated aneuploidy in development, ageing and cancer. *Nat Rev Genet* 2024. PMID 39169218. doi:10.1038/s41576-024-00762-6.
22. García-Castillo H et al. Clinical and genetic heterogeneity in patients with mosaic variegated aneuploidy: delineation of clinical subtypes. *Am J Med Genet A* 2008. PMID 18548531. doi:10.1002/ajmg.a.32315.
23. Limwongse C et al. Child with mosaic variegated aneuploidy and embryonal rhabdomyosarcoma. *Am J Med Genet* 1999. PMID 9916837.
24. McLaren W et al. The Ensembl Variant Effect Predictor. *Genome Biol* 2016. PMID 27268795. doi:10.1186/s13059-016-0974-4.
25. Cheng J et al. Accurate proteome-wide missense variant effect prediction with AlphaMissense. *Science* 2023. PMID 37733863. doi:10.1126/science.adg7492.
26. Ioannidis NM et al. REVEL: An Ensemble Method for Predicting the Pathogenicity of Rare Missense Variants. *Am J Hum Genet* 2016. PMID 27666373. doi:10.1016/j.ajhg.2016.08.016.
27. Rentzsch P et al. CADD: predicting the deleteriousness of variants throughout the human genome. *Nucleic Acids Res* 2019. PMID 30371827. doi:10.1093/nar/gky1016.
28. Jaganathan K et al. Predicting Splicing from Primary Sequence with Deep Learning. *Cell* 2019. PMID 30661751. doi:10.1016/j.cell.2018.12.015.
29. Lindeboom RG et al. The rules and impact of nonsense-mediated mRNA decay in human cancers. *Nat Genet* 2016. PMID 27618451. doi:10.1038/ng.3664.
30. Richards S et al. Standards and guidelines for the interpretation of sequence variants: a joint consensus recommendation of the American College of Medical Genetics and Genomics and the Association for Molecular Pathology. *Genet Med* 2015. PMID 25741868. doi:10.1038/gim.2015.30.
31. Jumper J et al. Highly accurate protein structure prediction with AlphaFold. *Nature* 2021. PMID 34265844. doi:10.1038/s41586-021-03819-2.
32. Varadi M et al. AlphaFold Protein Structure Database in 2024: providing structure coverage for over 214 million protein sequences. *Nucleic Acids Res* 2024. PMID 37933859. doi:10.1093/nar/gkad1011.
33. Dieckhaus H et al. Transfer learning to leverage larger datasets for improved prediction of protein stability changes. *Proc Natl Acad Sci U S A* 2024. PMID 38285937. doi:10.1073/pnas.2314853121.
34. Lin Z et al. Evolutionary-scale prediction of atomic-level protein structure with a language model. *Science* 2023. PMID 36927031. doi:10.1126/science.ade2574.
35. Rives A et al. Biological structure and function emerge from scaling unsupervised learning to 250 million protein sequences. *Proc Natl Acad Sci U S A* 2021. PMID 33876751. doi:10.1073/pnas.2016239118.
36. McWilliam SJ et al. Aminoglycoside-induced nephrotoxicity in children. *Pediatr Nephrol* 2017. PMID 27848094. doi:10.1007/s00467-016-3533-z.
37. Nishitani-Isa M et al. Rhabdomyosarcoma with premature chromatid separation-mosaic variegated aneuploidy syndrome: Reduced-intensity chemotherapy. *Pediatr Int* 2019. PMID 31184400. doi:10.1111/ped.13849.
38. Kirkegaard T et al. Heat shock protein-based therapy as a potential candidate for treating the sphingolipidoses. *Sci Transl Med* 2016. PMID 27605553. doi:10.1126/scitranslmed.aad9823.
39. Manuvakhova M et al. Aminoglycoside antibiotics mediate context-dependent suppression of termination codons in a mammalian translation system. *RNA* 2000. PMID 10917599. doi:10.1017/s1355838200000716.
40. Dai W et al. Slippage of mitotic arrest and enhanced tumor development in mice with BubR1 haploinsufficiency. *Cancer Res* 2004. PMID 14744753. doi:10.1158/0008-5472.can-03-3119.
41. Polier S et al. ATP-competitive inhibitors block protein kinase recruitment to the Hsp90-Cdc37 system. *Nat Chem Biol* 2013. PMID 23502424. doi:10.1038/nchembio.1212.
42. Moriyama B et al. Adverse interactions between antifungal azoles and vincristine: review and analysis of cases. *Mycoses* 2012. PMID 22126626. doi:10.1111/j.1439-0507.2011.02158.x.
43. Laberge RM et al. MTOR regulates the pro-tumorigenic senescence-associated secretory phenotype by promoting IL1A translation. *Nat Cell Biol* 2015. PMID 26147250. doi:10.1038/ncb3195.
44. Goldbach-Mansky R et al. Neonatal-onset multisystem inflammatory disease responsive to interleukin-1beta inhibition. *N Engl J Med* 2006. PMID 16899778. doi:10.1056/NEJMoa055137.
45. Sibley CH et al. Sustained response and prevention of damage progression in patients with neonatal-onset multisystem inflammatory disease treated with anakinra: a cohort study to determine three- and five-year outcomes. *Arthritis Rheum* 2012. PMID 22294344. doi:10.1002/art.34409.
46. Bali P et al. Inhibition of histone deacetylase 6 acetylates and disrupts the chaperone function of heat shock protein 90: a novel basis for antileukemia activity of histone deacetylase inhibitors. *J Biol Chem* 2005. PMID 15937340. doi:10.1074/jbc.C500186200.
47. Pirinen E et al. Niacin Cures Systemic NAD(+) Deficiency and Improves Muscle Performance in Adult-Onset Mitochondrial Myopathy. *Cell Metab* 2020. PMID 32386566. doi:10.1016/j.cmet.2020.04.008.
48. Aaltio J et al. Niacin supplementation in a child with novel MTTN variant m.5670A>G causing early onset mitochondrial myopathy and NAD(+) deficiency. *Neuromuscul Disord* 2024. PMID 39173541. doi:10.1016/j.nmd.2024.07.005.
49. Avalos JL et al. Mechanism of sirtuin inhibition by nicotinamide: altering the NAD(+) cosubstrate specificity of a Sir2 enzyme. *Mol Cell* 2005. PMID 15780941. doi:10.1016/j.molcel.2005.02.022.
50. Wu LE et al. SIRT2 transgenic over-expression does not impact lifespan in mice. *Aging Cell* 2023. PMID 38009412. doi:10.1111/acel.14027.
51. Colletti RB et al. Niacin treatment of hypercholesterolemia in children. *Pediatrics* 1993. PMID 8516088.
52. Mengel E et al. Efficacy and safety of arimoclomol in Niemann-Pick disease type C: Results from a double-blind, randomised, placebo-controlled, multinational phase 2/3 trial of a novel treatment. *J Inherit Metab Dis* 2021. PMID 34418116. doi:10.1002/jimd.12428.
53. Gray J et al. Heat shock protein amplification improves cerebellar myelination in the Npc1(nih) mouse model. *EBioMedicine* 2022. PMID 36455410. doi:10.1016/j.ebiom.2022.104374.
54. Tveten K et al. 4-Phenylbutyrate restores the functionality of a misfolded mutant low-density lipoprotein receptor. *FEBS J* 2007. PMID 17408384. doi:10.1111/j.1742-4658.2007.05735.x.
55. Choi E et al. BubR1 acetylation at prometaphase is required for modulating APC/C activity and timing of mitosis. *EMBO J* 2009. PMID 19407811. doi:10.1038/emboj.2009.123.
56. Mannick JB et al. mTOR inhibition improves immune function in the elderly. *Sci Transl Med* 2014. PMID 25540326. doi:10.1126/scitranslmed.3009892.
57. Mannick JB et al. TORC1 inhibition enhances immune function and reduces infections in the elderly. *Sci Transl Med* 2018. PMID 29997249. doi:10.1126/scitranslmed.aaq1564.
58. Mascarenhas L et al. Randomized Phase II Trial of Bevacizumab or Temsirolimus in Combination With Chemotherapy for First Relapse Rhabdomyosarcoma: A Report From the Children's Oncology Group. *J Clin Oncol* 2019. PMID 31513481. doi:10.1200/JCO.19.00576.
59. Gupta AA et al. Addition of temsirolimus to chemotherapy in children, adolescents, and young adults with intermediate-risk rhabdomyosarcoma (ARST1431): a randomised, open-label, phase 3 trial from the Children's Oncology Group. *Lancet Oncol* 2024. PMID 38936378. doi:10.1016/S1470-2045(24)00255-9.
60. Zhu Y et al. The Achilles' heel of senescent cells: from transcriptome to senolytic drugs. *Aging Cell* 2015. PMID 25754370. doi:10.1111/acel.12344.
61. Hickson LJ et al. Senolytics decrease senescent cells in humans: Preliminary report from a clinical trial of Dasatinib plus Quercetin in individuals with diabetic kidney disease. *EBioMedicine* 2019. PMID 31542391. doi:10.1016/j.ebiom.2019.08.069.
62. Gore L et al. Dasatinib in Pediatric Patients With Chronic Myeloid Leukemia in Chronic Phase: Results From a Phase II Trial. *J Clin Oncol* 2018. PMID 29498925. doi:10.1200/JCO.2017.75.9597.
63. Yousefzadeh MJ et al. Fisetin is a senotherapeutic that extends health and lifespan. *EBioMedicine* 2018. PMID 30279143. doi:10.1016/j.ebiom.2018.09.015.
64. Tang YC et al. Identification of aneuploidy-selective antiproliferation compounds. *Cell* 2011. PMID 21315436. doi:10.1016/j.cell.2011.01.017.
65. Mosallaei D et al. Molecular and Clinical Outcomes After Intravenous Gentamicin Treatment for Patients With Junctional Epidermolysis Bullosa Caused by Nonsense Variants. *JAMA Dermatol* 2022. PMID 35234826. doi:10.1001/jamadermatol.2021.5992.
66. Xu M et al. JAK inhibition alleviates the cellular senescence-associated secretory phenotype and frailty in old age. *Proc Natl Acad Sci U S A* 2015. PMID 26578790. doi:10.1073/pnas.1515386112.
67. Ochiai K et al. Long-term remission of bilateral Wilms tumors that developed from premature separation of chromatids/mosaic variegated aneuploidy syndrome due to bilateral nephrectomy and peritoneal dialysis. *Pediatr Blood Cancer* 2019. PMID 31081598. doi:10.1002/pbc.27804.
68. Hoy SM et al. Pimitespib: First Approval. *Drugs* 2022. PMID 35986838. doi:10.1007/s40265-022-01764-6.
