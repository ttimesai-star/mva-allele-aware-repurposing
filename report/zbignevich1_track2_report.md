# Allele-aware drug repurposing for BUB1B mosaic variegated aneuploidy (MVA1)

*Rare Disease, Real Kid: MVA Hackathon 2026 - Track 2 (Drug Repurposing) - submission by zbignevich1*


> **Nature of this document.** This report proposes *research hypotheses requiring clinical evaluation*. It is not
> medical advice and not a treatment recommendation for any individual. It contains no doses, schedules or
> instructions for use in children. Every candidate below is an approved medicine proposed for laboratory and
> clinical evaluation by qualified investigators and the treating team; several carry serious toxicities.
> The controlled hackathon dataset was **not** used: the framework is built from public data only and is
> parametrised by genotype so that the panel or treating team can instantiate it locally with the confirmed variants.

**AI disclosure (required).** Anthropic Claude Opus 5.5 via Claude Code; [AI DISCLOSURE: tier/opt-out - to confirm].
No hackathon (controlled) data was ever provided to the AI system.

## Summary

Mosaic variegated aneuploidy type 1 (MVA1) is caused by biallelic variants in *BUB1B*, which encodes the spindle
assembly checkpoint kinase BUBR1 [1]. We argue that the right unit for drug repurposing in MVA1 is not the
gene but the **allele class**. In patient cells, MVA1 mutations fall into two functional classes: those that lower
BUBR1 abundance (truncating alleles give no transcript; missense alleles in or near the kinase domain make protein
that is turned over fast) and those with specific checkpoint or attachment defects [2]. A
dosage-raising drug can only help if at least one allele still makes protein. We therefore built a transparent,
citation-locked scoring framework that (i) classifies each allele, (ii) maps four mechanism-anchored intervention
strategies onto the disease chain, (iii) scores ten approved medicines on six evidence dimensions, (iv) ranks them
separately for five biallelic genotype configurations and (v) tests rank robustness over 5,000 perturbed weightings.

For the most common MVA1 configuration (one truncating plus one kinase-domain missense allele), the leading
hypotheses are: **H1** NAD+ repletion with niacin to support SIRT2-dependent BUBR1 stability (rank 1 in 88% and top
three in 98% of perturbed rankings); **H2** heat-shock-response amplification with arimoclomol to stabilise a
folding-defective missense BUBR1; **H3** low-dose intermittent mTORC1 dampening (everolimus) against progeroid muscle
decline; and **H4** intermittent senolytic clearance (dasatinib, with quercetin as a non-licensed partner) of
senescent, often aneuploid, cells. The ranking changes with genotype: if both alleles are null or the missense allele
is stable, H1 and H2 drop and H3/H4 lead. A machine-generated **do-no-harm layer** lists approved medicines whose
mechanism depends on an intact spindle checkpoint (23 approved tubulin-targeting agents in ChEMBL) and other
clinical cautions for the treating team. All 47 references were verified against PubMed (title and DOI) by the
pipeline itself; the build fails if any citation is missing or unverifiable.

## 1. Variant mechanism: what goes wrong, and where a medicine could act

**Loss of function by dosage.** MVA1 is recessive and biallelic; the shared consequence of the two alleles is too
little functional BUBR1 [1, 2]. In biallelic patient cell lines, a missense allele usually
pairs with a truncating allele, and the cells show an impaired mitotic checkpoint, chromosome alignment defects and
low overall BUBR1; re-expressing BUBR1 restores checkpoint activity, which establishes BUBR1 dysfunction as causal
[2]. We characterise the disease as **loss of function with a dosage component**; in the sources we
reviewed, no evidence supports a gain-of-function mechanism for MVA1 alleles.

**Two allele classes, two kinds of missense.** Truncating alleles produce no transcript; missense alleles in or near
the kinase domain produce protein with high turnover; other missense alleles keep normal abundance but impair
checkpoint activity or microtubule attachment [2]. Selected missense mutants (I909T, L1012P) depend
on HSP90 to escape proteasomal degradation, i.e. they behave as folding-defective proteins [2].
Regulatory alleles that lower expression of wild-type protein have also been described, including an intergenic
variant upstream of *BUB1B* in Japanese families [3, 4].

**Pathway and downstream consequences** (Figure 1; evidence level in brackets: 3 = shown in BubR1-insufficient mice or
MVA patient cells, 2 = shown in cells or one model genotype, 1 = inference):

1. Low BUBR1 -> weak spindle checkpoint and attachment errors [3] [2, 5].
2. -> chromosome missegregation and mosaic aneuploidy [3] [1, 6].
3. -> p53-p21 and p16Ink4a senescence with a senescence-associated secretory phenotype (SASP) [3]
   [7, 8, 9].
4. -> progeroid tissue decline: cataract, lipodystrophy, sarcopenia [3] [6, 10].
5. -> embryonal cancers (rhabdomyosarcoma, Wilms tumour, leukaemia) [3] [1, 11].
6. Allele-specific branch: heterozygous mice carrying the truncating X753 allele show skeletal-muscle mTORC1
   hyperactivity (raised p70 S6 kinase and 4EBP1 phosphorylation) that tracks with sarcopenia, while carriers of a
   missense allele with a similar BUBR1 reduction do not [2] [12].
7. Counter-regulation: NAD+-dependent SIRT2 deacetylation maintains BUBR1 levels; NAD+ decline lowers BUBR1, and the
   NAD+ precursor NMN raised BUBR1 abundance in vivo [3] [13, 14].
8. Aneuploid, arrested cells can be cleared by NK cells via an NF-kB-driven program [2] [9, 15].

The genetics are broader than *BUB1B*: biallelic *CEP57*, *TRIP13*, *BUB1* and *MAD1L1* variants cause related
aneuploidy syndromes [16, 17, 18, 19, 20]. Clinically, MVA
is heterogeneous, with distinct subtypes and cancer predisposition [21, 22].
Increased BubR1 expression protects mice against aneuploidy and cancer and extends healthy lifespan [5],
which makes **raising residual BUBR1** the most causal lever and **limiting downstream senescence/SASP** the most
allele-independent one.

![Figure 1. Mechanism map: the MVA1 disease chain and the intervention points of strategies S1-S3 (S4 acts on aneuploid clones). Solid arrows: level-3 evidence; dashed: one mouse genotype. The do-no-harm layer is shown at the bottom.](https://github.com/ttimesai-star/mva-allele-aware-repurposing/raw/main/figures/mechanism_map.png)

## 2. Approach

The pipeline (`run_all.py`) is fully scripted and reruns offline from a committed HTTP cache.

1. **Reference verification** (`verify_refs.py`). Each of 47 references is resolved against PubMed E-utilities; title
   word overlap and DOI must match. Any citation key used in the knowledge base but not defined, or any reference that
   fails verification, stops the build. Result: 47/47 verified (`outputs/reference_verification.csv`).
2. **Allele-class landscape** (`clinvar_landscape.py`). All single-gene Pathogenic/Likely pathogenic ClinVar records for
   seven MVA genes are classified as truncating (frameshift/splice), nonsense, kinase-window missense (unstable),
   other missense (functional) or regulatory, using the UniProt O60566 kinase domain (766-1050) extended N-terminally to
   residue 720 so that R727C, which is unstable in patient cells, counts as "near" [2].
3. **Rule check against measured protein abundance.** The kinase-window rule is compared with the eight MVA1 missense
   alleles whose protein abundance was measured [2].
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
8. **Genotype instantiation** (`classify_variants.py`). Given the two confirmed alleles, the script prints their classes
   and the matching ranked list locally; the genotype never needs to enter the repository.

## 3. Results

### 3.1 The dosage-responsive allele class is invisible in ClinVar

Of 78 single-gene Pathogenic/Likely pathogenic *BUB1B* records, 76 are truncating (40 frameshift/splice, 36 nonsense),
one is a missense outside the kinase window and one is unclassified (`outputs/clinvar_allele_class_summary.csv`). In
contrast, of 421 *BUB1B* missense records linked to MVA or premature chromatid separation, 392 are of uncertain
significance and 22 have conflicting classifications; one is Likely pathogenic
(`outputs/bub1b_missense_clinvar_classes.csv`). The missense half of a typical MVA1 genotype - exactly the half that
determines whether a dosage-raising drug can work - is therefore usually a variant of uncertain significance in public
databases. The kinase-window rule reproduces measured protein abundance for 7 of 8 published MVA1 missense alleles;
the single discordant allele (Q921H, normal abundance) occurs in cis with the unstable L844F allele
[2] (`outputs/kinase_window_rule_check.csv`, Figure 2).

![Figure 2. BUBR1 allele map. Above: ClinVar Pathogenic/Likely pathogenic variants with a protein position. Below: published MVA1 missense alleles coloured by measured protein abundance; the shaded window marks the "in or near kinase domain" rule.](https://github.com/ttimesai-star/mva-allele-aware-repurposing/raw/main/figures/bubr1_alleles.png)

### 3.2 The ranking depends on the genotype

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
| gentamicin (S1c) | 0.10 | 0.42 | 0.10 | 0.10 | 0.10 |

Table 1. Allele-aware scores (`outputs/ranking_by_genotype.csv`; Figure 3). The first column is the most common MVA1
configuration [2] and our primary scenario. Scores rank hypotheses for testing; they are not
predictions of efficacy.

![Figure 3. Allele-aware score by BUB1B genotype configuration.](https://github.com/ttimesai-star/mva-allele-aware-repurposing/raw/main/figures/ranking_heatmap.png)

### 3.3 Robustness

Under 5,000 perturbations of weights and expert scores, niacin ranked first in 88.2% and in the top three in 97.7% of
draws for the primary genotype. Arimoclomol (54.7% top three) and everolimus (52.9%) compete for second place,
followed by dasatinib (33.7%); ruxolitinib and gentamicin never entered the top three
(`outputs/rank_stability_G1_TRUNC_x_MSUNSTABLE.csv`, Figure 4). The top of the list is robust to reasonable
disagreement about weights; places two to four are not, and we present them as a group.

![Figure 4. Rank stability for the primary genotype across 5,000 perturbed rankings.](https://github.com/ttimesai-star/mva-allele-aware-repurposing/raw/main/figures/rank_stability.png)

## 4. Repurposing hypotheses

Evidence tier for each hypothesis: **A** - the intervention itself was tested in a BubR1/MVA model; **B** - the
mechanism was validated (genetically or with a related agent) in a BubR1/MVA model and the drug engages that
mechanism in humans; **C** - mechanistic inference from other systems; **D** - exploratory. None of the hypotheses
below has tier A evidence.

### H1. NAD+ repletion with niacin to stabilise residual BUBR1 (tier B; strategy S1)

*Rationale.* BUBR1 declines with NAD+; SIRT2 keeps BUBR1 deacetylated and abundant; the NAD+ precursor NMN raised
BUBR1 abundance in vivo, and SIRT2 overexpression extended median lifespan of BubR1 hypomorphic mice [13];
SIRT2 deacetylation also governs BUBR1 degradation during mitosis [14]. More BUBR1 protects against
aneuploidy and cancer [5]. Niacin is a licensed NAD+ precursor that raised blood NAD+ up to 8-fold and
restored muscle NAD+ to control levels in adult patients with mitochondrial myopathy [23] and was used in a
child with NAD+ deficiency [24]. We prefer niacin to nicotinamide because nicotinamide is itself a sirtuin
inhibitor [25].

*Genotype dependence.* Only meaningful if at least one allele makes protein (unstable missense, regulatory or
functional missense); predicted to be inert for two null alleles.

*Counter-evidence and gaps.* Niacin itself has never been tested on BUBR1; SIRT2 overexpression did not extend lifespan
in wild-type mice [26]; paediatric niacin experience (hypercholesterolaemia) showed frequent reversible adverse
effects, including aminotransferase elevations [27], and the US label states that paediatric safety is not
established. Raising NAD+ in a cancer-predisposed child has theoretical tumour-metabolism implications.

*Decisive test.* BUBR1 western blot, nocodazole checkpoint assay and mitotic-error rate in MVA1 patient fibroblasts or
lymphoblastoid cells with and without an NAD+ precursor, stratified by allele class.

### H2. Heat-shock-response amplification with arimoclomol to rescue a folding-defective missense BUBR1 (tier C; S1b)

*Rationale.* Kinase-domain missense mutants turn over about twice as fast as wild type in patient cells and depend on
HSP90 to avoid proteasomal clearance [2]. Arimoclomol amplifies the heat-shock response to counter
protein misfolding; in a randomised placebo-controlled trial in patients aged 2-18 years with Niemann-Pick disease
type C it reduced annual disease progression by 65% [28], and heat-shock protein amplification improved
cerebellar myelination in an NPC mouse model [29]. Sodium phenylbutyrate, a licensed chemical chaperone, is a
weaker alternative: its rescue evidence is for a secretory protein (LDLR G544V) [30], whereas BUBR1 is
cytosolic, and it also inhibits histone deacetylases, which could alter BUBR1 acetylation [31].

*Genotype dependence.* Only for unstable missense alleles; near-zero fit otherwise (Table 1).

*Counter-evidence and gaps.* No data on BUBR1; heat-shock proteins also support tumour cells (theoretical);
treatment-related urticaria and angioedema were reported [28].

*Decisive test.* Cycloheximide chase of BUBR1 in patient cells carrying a kinase-domain missense allele, with and
without arimoclomol; readout BUBR1 half-life and checkpoint competence.

### H3. Low-dose intermittent mTORC1 dampening (everolimus; sirolimus as alternative) against progeroid decline (tier B-/C; S2)

*Rationale.* Truncating-allele (X753) carrier mice show skeletal-muscle mTORC1 hyperactivity that tracks with
sarcopenia [12]; mTOR drives the pro-tumorigenic SASP via IL1A translation [32]. Low-dose
everolimus regimens improved, rather than suppressed, influenza-vaccine responses and reduced infections in older
adults [33, 34]. Everolimus carries a paediatric label from age 1 (TSC-associated SEGA; openFDA).

*Scope.* Proposed for tissue ageing (muscle, fat), not as an anticancer drug: temsirolimus had activity in relapsed
rhabdomyosarcoma [35] but did not improve event-free survival when added to frontline therapy for
intermediate-risk rhabdomyosarcoma [36].

*Counter-evidence and gaps.* Rapamycin analogues have not been tested in BubR1 mice in the sources reviewed;
immunosuppression could blunt NK-cell clearance of aneuploid cells [15].

*Decisive test.* p-S6K/p-4EBP1 and SASP panel in MVA fibroblasts; muscle function in BubR1 +/X753 mice.

### H4. Intermittent senolytic clearance (dasatinib, with quercetin) of senescent aneuploid cells (tier B; S3)

*Rationale.* In BubR1 progeroid mice, genetic removal of p16Ink4a-positive cells delayed cataract, lipodystrophy and
sarcopenia, and late-life clearance attenuated established disease [10]; p16 and p21 shape these phenotypes
[7, 8]. Dasatinib plus quercetin is the senolytic pair identified by transcriptome-guided
screening [37] and reduced senescent-cell burden in a first human study [38]. Dasatinib has
paediatric experience in chronic myeloid leukaemia [39]. This strategy is allele-independent.

*Counter-evidence and gaps.* Senolytic drugs were not tested in BubR1 mice in our sources; dasatinib causes
myelosuppression and effusions; quercetin is a supplement, not a licensed medicine (fisetin is another non-licensed
senolytic [40]). Arrested aneuploid cells may already be cleared by the immune system
[9, 15], which lowers the expected added benefit.

*Decisive test.* SA-beta-gal and p16 counts in MVA fibroblast cultures after pulsed exposure; cataract and fat
endpoints in BubR1 H/H mice.

### H5. Exploratory hypotheses (tier C-D)

* **Aneuploid-cell stress (S4).** Chloroquine (autophagy inhibition) and AICAR (energy stress) were selectively
  anti-proliferative in trisomic mouse cells and aneuploid cancer lines [41]. Hydroxychloroquine and metformin
  (an AMPK activator; inference only) are licensed proxies. No MVA data.
* **Read-through of a nonsense allele (S1c).** Intravenous gentamicin produced read-through in children with nonsense
  variants in another disease [42]. Relevant only if one allele is a premature stop codon (36 of 78
  ClinVar P/LP *BUB1B* records); nephro- and ototoxicity make this a cell-assay probe, not a long-term option.
* **Ruxolitinib.** JAK inhibition suppressed the SASP and frailty in old mice [43]; weakest link to BUBR1 and a
  high infection risk; ranked last but one.

## 5. Do-no-harm layer: risks for the clinician to weigh

These flags are generated from mechanism and literature. They are research observations, not prescribing advice,
and are **not** a recommendation to withhold standard-of-care oncology treatment, which remains the decision of the
treating oncology team.

1. **Spindle-poison antimitotics.** BubR1-insufficient cells are defective in checkpoint activation and slip out of
   mitotic arrest [44]. The automated scan lists 23 approved tubulin-targeting agents in ChEMBL
   (`outputs/liability_scan.csv`), including vinca alkaloids and taxanes used in paediatric oncology, as well as
   benzimidazole anthelmintics (albendazole, mebendazole, thiabendazole) and griseofulvin, which bind tubulin.
   Whether the latter matter clinically in MVA is unknown; we flag them because they are commonly used in children.
2. **Severe treatment toxicity.** Case reports describe MVA1 children with Wilms tumour and with rhabdomyosarcoma, the
   latter managed with reduced-intensity chemotherapy [45, 46]; genomic instability disorders
   can show exquisite sensitivity to chemotherapy and ionising radiation [11].
3. **Radiation in surveillance.** The AACR workshop update supports renal ultrasound surveillance every 3 months from
   birth until age 7 for all MVA conditions and avoidance of radiation exposure [11].
4. **HSP90 inhibitors.** Kinase-domain missense BUBR1 needs HSP90 to escape degradation [2]; HSP90
   inhibition would be expected to lower it further. Pimitespib is approved in Japan [47]. An HSP90 inhibitor
   (17-AAG) was also among the aneuploidy-selective compounds [41], so this is a genuine trade-off for the
   clinician, not a settled contraindication.
5. **Mitotic kinase and Eg5 inhibitors** (Aurora, Polo-like kinase, KIF11): no approved agents matched in ChEMBL at the
   time of the scan; flagged for investigational use.
6. **Candidate-specific risks.** Immunosuppression with mTOR/JAK inhibitors (could blunt NK clearance of aneuploid
   cells [15]); myelosuppression with dasatinib; hepatotoxicity and flushing with niacin [27];
   nephro- and ototoxicity with gentamicin [42]; retinal toxicity with long-term hydroxychloroquine
   (label).
7. **Cancer predisposition.** Any agent that supports cell survival (NAD+ repletion, heat-shock proteins) must be
   weighed against tumour risk in a cancer-predisposition syndrome before any clinical use.

## 6. Limitations

* **The child's genotype was not used.** We did not download the controlled dataset. The Track 1 scoring rules imply
  a pair of variants; our framework covers five plausible biallelic configurations and can be instantiated locally
  (`classify_variants.py`). If both alleles are null, H1 and H2 do not apply.
* **Expert scores are judgements.** Every score is referenced and the sensitivity analysis shows which conclusions
  survive disagreement, but the scores were assigned by one participant with AI assistance and need independent review.
* **The kinase-window rule rests on eight measured alleles** from one study [2]; a functional assay
  (cycloheximide chase) should replace the rule for any new allele.
* **Literature bias.** The candidate list comes from mechanism-first literature review, not an unbiased screen; drugs
  acting through mechanisms outside the four strategies would be missed.
* **No efficacy data in MVA for any candidate.** Mouse evidence (BubR1 hypomorphs, X753 and L1002P alleles) may not
  translate; some benefits were sex-specific [13] or not reproduced in wild-type mice [26].
* **Database snapshot.** ClinVar, ChEMBL, openFDA and Open Targets were queried in October 2026; the HTTP cache makes
  the snapshot reproducible but not current.

## 7. Scalability

The framework separates generic machinery (reference verification, ClinVar allele-class landscape, ChEMBL/openFDA
annotation, liability scan, scoring, sensitivity analysis) from a disease-specific YAML knowledge base. It already
profiles six further MVA genes (*TRIP13*, *CEP57*, *BUB1*, *MAD1L1*, *CENATAC*, *BUB3*), and Open Targets lists 11
MVA-associated genes (`outputs/mva_gene_associations.csv`). Adding a disease means writing a new `evidence.yaml`;
adding a patient means two command-line arguments. The pipeline runs in minutes on a laptop without a GPU and reruns
offline.

## 8. Reproducibility

`python run_all.py` reproduces every table and figure (offline from the committed cache after the first run);
`pytest` covers allele classification, citation completeness, genotype fit and deterministic ranking. Code is
MIT-licensed; submissions are released under CC BY 4.0 per the hackathon rules. Repository:
https://github.com/ttimesai-star/mva-allele-aware-repurposing (private until the hackathon closes, then public).

## 9. Methods description (Track 2 template)

**Team name.** zbignevich1 (individual participant).

**Approach in detail.** Mechanism-first, allele-aware repurposing: (1) characterise the BUB1B disease chain and the
two functional allele classes from primary literature; (2) define four intervention strategies anchored to nodes of
that chain (raise residual BUBR1; help unstable missense BUBR1 fold; dampen mTORC1/SASP; clear senescent cells) plus
exploratory ones; (3) select approved medicines acting on each node and confirm approval status (ChEMBL max phase 4)
and paediatric label text (openFDA); (4) score each candidate on six referenced evidence dimensions; (5) weight by
allele fit for five genotype configurations; (6) test rank stability by Monte Carlo perturbation; (7) scan all
approved drugs for checkpoint-dependent mechanisms to build the do-no-harm layer.

**Generative AI (required).** Anthropic Claude Opus 5.5 via Claude Code; [AI DISCLOSURE: tier/opt-out - to confirm].
The AI assisted with literature-search queries, code and drafting; all references were verified programmatically
against PubMed and key claims were checked against abstracts or full texts. No controlled hackathon data was given to
the AI.

**Automated or manual identification.** Hybrid. Candidate selection and expert scores come from manual, AI-assisted
literature review and curation; database annotation, reference verification, the ClinVar landscape, the liability
scan, scoring and sensitivity analysis are automated.

**Manual literature review and curation.** Targeted PubMed searches on BUBR1 biology, BubR1 mouse models, NAD+/SIRT2,
senescence and senolytics, mTOR and SASP, aneuploidy-selective compounds, proteostasis drugs and the paediatric safety
of each candidate. Each mechanism edge, strategy and score in `config/evidence.yaml` carries citation keys;
counter-evidence is recorded explicitly (`counter_refs`).

**Public data only?** Yes. Only public databases and literature were used.

**Public data sources.** PubMed/E-utilities (literature and reference verification), ClinVar (variant classes),
UniProt O60566 (domain map), ChEMBL (approval status, mechanisms, liability scan), openFDA drug labels (paediatric use),
Open Targets Platform (MVA-associated genes), MONDO (MONDO_0009759).

**Proprietary data.** None.

**Mechanism characterisation.** Loss of function with a dosage component (Section 1): truncating alleles give no
transcript, kinase-domain missense alleles give unstable protein, other missense alleles give specific checkpoint or
attachment defects [2]; pathway: spindle assembly checkpoint -> chromosome missegregation ->
p53/p21/p16 senescence and SASP -> progeroid decline and embryonal cancer, with an allele-specific mTORC1 branch.

**Time and effort.** About 20 person-hours including AI-assisted work; pipeline runtime of minutes on a laptop; no
paid compute.

**Method abstract.** See *Summary*. Strengths: allele-aware rather than gene-level ranking; every claim tied to a
verified reference; explicit counter-evidence; quantified robustness; machine-generated safety layer; public data only;
scales by configuration. Limitations: Section 6.

## Acknowledgement

This work was made possible through the Hackathon, organized by Sage Bionetworks in partnership with the MVA Society,
Hugging Face, and BEACON (The Benchmarking, Evaluation, and Assessment Consortium for Science), with prize sponsorship
from AWS and Anthropic. We are deeply grateful to the child and their family who generously contributed their data and
their story to advance research into this rare disease. We acknowledge their trust in making this Hackathon possible.

## References

1. Hanks S et al. Constitutional aneuploidy and cancer predisposition caused by biallelic mutations in BUB1B. *Nat Genet* 2004. PMID 15475955. doi:10.1038/ng1449.
2. Suijkerbuijk SJ et al. Molecular causes for BUBR1 dysfunction in the human cancer predisposition syndrome mosaic variegated aneuploidy. *Cancer Res* 2010. PMID 20516114. doi:10.1158/0008-5472.CAN-09-4319.
3. Ochiai H et al. TALEN-mediated single-base-pair editing identification of an intergenic mutation upstream of BUB1B as causative of PCS (MVA) syndrome. *Proc Natl Acad Sci U S A* 2014. PMID 24344301. doi:10.1073/pnas.1317008111.
4. Matsuura S et al. Monoallelic BUB1B mutations and defective mitotic-spindle checkpoint in seven families with premature chromatid separation (PCS) syndrome. *Am J Med Genet A* 2006. PMID 16411201. doi:10.1002/ajmg.a.31069.
5. Baker DJ et al. Increased expression of BubR1 protects against aneuploidy and cancer and extends healthy lifespan. *Nat Cell Biol* 2013. PMID 23242215. doi:10.1038/ncb2643.
6. Baker DJ et al. BubR1 insufficiency causes early onset of aging-associated phenotypes and infertility in mice. *Nat Genet* 2004. PMID 15208629. doi:10.1038/ng1382.
7. Baker DJ et al. Opposing roles for p16Ink4a and p19Arf in senescence and ageing caused by BubR1 insufficiency. *Nat Cell Biol* 2008. PMID 18516091. doi:10.1038/ncb1744.
8. Baker DJ et al. p21 both attenuates and drives senescence and aging in BubR1 progeroid mice. *Cell Rep* 2013. PMID 23602569. doi:10.1016/j.celrep.2013.03.028.
9. Santaguida S et al. Chromosome Mis-segregation Generates Cell-Cycle-Arrested Cells with Complex Karyotypes that Are Eliminated by the Immune System. *Dev Cell* 2017. PMID 28633018. doi:10.1016/j.devcel.2017.05.022.
10. Baker DJ et al. Clearance of p16Ink4a-positive senescent cells delays ageing-associated disorders. *Nature* 2011. PMID 22048312. doi:10.1038/nature10600.
11. Nakano Y et al. Update on Recommendations for Cancer Screening and Surveillance in Children with Genomic Instability Disorders. *Clin Cancer Res* 2024. PMID 39264246. doi:10.1158/1078-0432.CCR-24-1098.
12. Sieben CJ et al. BubR1 allelic effects drive phenotypic heterogeneity in mosaic-variegated aneuploidy progeria syndrome. *J Clin Invest* 2020. PMID 31738183. doi:10.1172/JCI126863.
13. North BJ et al. SIRT2 induces the checkpoint kinase BubR1 to increase lifespan. *EMBO J* 2014. PMID 24825348. doi:10.15252/embj.201386907.
14. Suematsu T et al. Deacetylation of the mitotic checkpoint protein BubR1 at lysine 250 by SIRT2 and subsequent effects on BubR1 degradation during the prometaphase/anaphase transition. *Biochem Biophys Res Commun* 2014. PMID 25285631. doi:10.1016/j.bbrc.2014.09.128.
15. Wang RW et al. Aneuploid senescent cells activate NF-κB to promote their immune clearance by NK cells. *EMBO Rep* 2021. PMID 34105235. doi:10.15252/embr.202052032.
16. Aziz K et al. Mosaic-variegated aneuploidy syndrome mutation or haploinsufficiency in Cep57 impairs tumor suppression. *J Clin Invest* 2018. PMID 30035751. doi:10.1172/JCI120316.
17. Yost S et al. Biallelic TRIP13 mutations predispose to Wilms tumor and chromosome missegregation. *Nat Genet* 2017. PMID 28553959. doi:10.1038/ng.3883.
18. Carvalhal S et al. Biallelic BUB1 mutations cause microcephaly, developmental delay, and variable effects on cohesion and chromosome segregation. *Sci Adv* 2022. PMID 35044816. doi:10.1126/sciadv.abk0114.
19. Villarroya-Beltri C et al. Biallelic germline mutations in MAD1L1 induce a syndrome of aneuploidy with high tumor susceptibility. *Sci Adv* 2022. PMID 36322655. doi:10.1126/sciadv.abq5914.
20. Malumbres M et al. Mosaic variegated aneuploidy in development, ageing and cancer. *Nat Rev Genet* 2024. PMID 39169218. doi:10.1038/s41576-024-00762-6.
21. García-Castillo H et al. Clinical and genetic heterogeneity in patients with mosaic variegated aneuploidy: delineation of clinical subtypes. *Am J Med Genet A* 2008. PMID 18548531. doi:10.1002/ajmg.a.32315.
22. Limwongse C et al. Child with mosaic variegated aneuploidy and embryonal rhabdomyosarcoma. *Am J Med Genet* 1999. PMID 9916837.
23. Pirinen E et al. Niacin Cures Systemic NAD(+) Deficiency and Improves Muscle Performance in Adult-Onset Mitochondrial Myopathy. *Cell Metab* 2020. PMID 32386566. doi:10.1016/j.cmet.2020.04.008.
24. Aaltio J et al. Niacin supplementation in a child with novel MTTN variant m.5670A>G causing early onset mitochondrial myopathy and NAD(+) deficiency. *Neuromuscul Disord* 2024. PMID 39173541. doi:10.1016/j.nmd.2024.07.005.
25. Avalos JL et al. Mechanism of sirtuin inhibition by nicotinamide: altering the NAD(+) cosubstrate specificity of a Sir2 enzyme. *Mol Cell* 2005. PMID 15780941. doi:10.1016/j.molcel.2005.02.022.
26. Wu LE et al. SIRT2 transgenic over-expression does not impact lifespan in mice. *Aging Cell* 2023. PMID 38009412. doi:10.1111/acel.14027.
27. Colletti RB et al. Niacin treatment of hypercholesterolemia in children. *Pediatrics* 1993. PMID 8516088.
28. Mengel E et al. Efficacy and safety of arimoclomol in Niemann-Pick disease type C: Results from a double-blind, randomised, placebo-controlled, multinational phase 2/3 trial of a novel treatment. *J Inherit Metab Dis* 2021. PMID 34418116. doi:10.1002/jimd.12428.
29. Gray J et al. Heat shock protein amplification improves cerebellar myelination in the Npc1(nih) mouse model. *EBioMedicine* 2022. PMID 36455410. doi:10.1016/j.ebiom.2022.104374.
30. Tveten K et al. 4-Phenylbutyrate restores the functionality of a misfolded mutant low-density lipoprotein receptor. *FEBS J* 2007. PMID 17408384. doi:10.1111/j.1742-4658.2007.05735.x.
31. Choi E et al. BubR1 acetylation at prometaphase is required for modulating APC/C activity and timing of mitosis. *EMBO J* 2009. PMID 19407811. doi:10.1038/emboj.2009.123.
32. Laberge RM et al. MTOR regulates the pro-tumorigenic senescence-associated secretory phenotype by promoting IL1A translation. *Nat Cell Biol* 2015. PMID 26147250. doi:10.1038/ncb3195.
33. Mannick JB et al. mTOR inhibition improves immune function in the elderly. *Sci Transl Med* 2014. PMID 25540326. doi:10.1126/scitranslmed.3009892.
34. Mannick JB et al. TORC1 inhibition enhances immune function and reduces infections in the elderly. *Sci Transl Med* 2018. PMID 29997249. doi:10.1126/scitranslmed.aaq1564.
35. Mascarenhas L et al. Randomized Phase II Trial of Bevacizumab or Temsirolimus in Combination With Chemotherapy for First Relapse Rhabdomyosarcoma: A Report From the Children's Oncology Group. *J Clin Oncol* 2019. PMID 31513481. doi:10.1200/JCO.19.00576.
36. Gupta AA et al. Addition of temsirolimus to chemotherapy in children, adolescents, and young adults with intermediate-risk rhabdomyosarcoma (ARST1431): a randomised, open-label, phase 3 trial from the Children's Oncology Group. *Lancet Oncol* 2024. PMID 38936378. doi:10.1016/S1470-2045(24)00255-9.
37. Zhu Y et al. The Achilles' heel of senescent cells: from transcriptome to senolytic drugs. *Aging Cell* 2015. PMID 25754370. doi:10.1111/acel.12344.
38. Hickson LJ et al. Senolytics decrease senescent cells in humans: Preliminary report from a clinical trial of Dasatinib plus Quercetin in individuals with diabetic kidney disease. *EBioMedicine* 2019. PMID 31542391. doi:10.1016/j.ebiom.2019.08.069.
39. Gore L et al. Dasatinib in Pediatric Patients With Chronic Myeloid Leukemia in Chronic Phase: Results From a Phase II Trial. *J Clin Oncol* 2018. PMID 29498925. doi:10.1200/JCO.2017.75.9597.
40. Yousefzadeh MJ et al. Fisetin is a senotherapeutic that extends health and lifespan. *EBioMedicine* 2018. PMID 30279143. doi:10.1016/j.ebiom.2018.09.015.
41. Tang YC et al. Identification of aneuploidy-selective antiproliferation compounds. *Cell* 2011. PMID 21315436. doi:10.1016/j.cell.2011.01.017.
42. Mosallaei D et al. Molecular and Clinical Outcomes After Intravenous Gentamicin Treatment for Patients With Junctional Epidermolysis Bullosa Caused by Nonsense Variants. *JAMA Dermatol* 2022. PMID 35234826. doi:10.1001/jamadermatol.2021.5992.
43. Xu M et al. JAK inhibition alleviates the cellular senescence-associated secretory phenotype and frailty in old age. *Proc Natl Acad Sci U S A* 2015. PMID 26578790. doi:10.1073/pnas.1515386112.
44. Dai W et al. Slippage of mitotic arrest and enhanced tumor development in mice with BubR1 haploinsufficiency. *Cancer Res* 2004. PMID 14744753. doi:10.1158/0008-5472.can-03-3119.
45. Ochiai K et al. Long-term remission of bilateral Wilms tumors that developed from premature separation of chromatids/mosaic variegated aneuploidy syndrome due to bilateral nephrectomy and peritoneal dialysis. *Pediatr Blood Cancer* 2019. PMID 31081598. doi:10.1002/pbc.27804.
46. Nishitani-Isa M et al. Rhabdomyosarcoma with premature chromatid separation-mosaic variegated aneuploidy syndrome: Reduced-intensity chemotherapy. *Pediatr Int* 2019. PMID 31184400. doi:10.1111/ped.13849.
47. Hoy SM et al. Pimitespib: First Approval. *Drugs* 2022. PMID 35986838. doi:10.1007/s40265-022-01764-6.
