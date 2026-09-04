# In vivo endogenous proximity proteomes of α-synuclein, Vps35 and Lrrk2 converge on the presynaptic vesicle module and identify Stxbp1 as a selective modifier of α-synuclein-dependent dopamine neuron loss

Analysis code and processed data for the computational figures, supplemental figures and
supplemental tables of the manuscript.

## Layout

```
data/                 processed, analysis-ready inputs, one copy of each file
  ibioid/             HiUGE-iBioID contrasts, proximity proteome sets, GO exports, replicates
  ibiocofrac/         co-fractionation tidy tables, modules, fractions, permutation results
  craave_seq/         screen library, counts, MAGeCK output, histology, western blot, MEA
  cortex/             cortical cohort contrasts
  human/              published human single-cell derivatives
  external/           third-party reference files
  PROVENANCE*.csv     md5, row count, source and consumers for every file above
src/                  modules imported by the builders
figures/main_01 .. main_07/       scripts/ and output/
figures/suppl_01 .. suppl_10/     scripts/ and output/
tables/               scripts/ and output/
MANIFEST.csv          one row per manuscript panel: deliverable, panel, kind, script, inputs, outputs
CHECKSUMS.csv         every produced file with its size and md5
```

## Environment

```
conda env create -f environment.yml && conda activate proximity-proteomics
```

or, into an existing Python 3.11 environment:

```
pip install -r requirements.txt
```

No R is required.

## Running

Every builder is standalone and writes into its own figure's `output/`:

```
python3 figures/main_02/scripts/venn_wt_shared.py
```

Three chains must be run in order, because a later step reads an earlier step's table:

```
figures/main_05/scripts/module_group_permutation.py    then module_group_fc_heatmap.py
                                                       and figures/suppl_06/scripts/module_fc_by_group.py
figures/main_06/scripts/module_profile_prepare.py      then module_permutation_bootstrap.py
                                                       then m11_fraction_profile.py
figures/suppl_06/scripts/module_coherence_permutation.py  then module_coherence_bh.py
                                                          then module_coherence_panel.py
figures/suppl_07/scripts/cortex_group_permutation.py   then panels_BCD.py
```

To run everything in that order:

```
python3 run_all.py
```

## Which script builds which panel

`MANIFEST.csv` has one row per panel of the manuscript, with `kind`:

- `script` - built by the listed script
- `manual` - immunoblot, microscopy or schematic; no computational producer
- `cytoscape` - network laid out in Cytoscape and finished in Illustrator; the node and edge
  tables are in `data/`, the layout is not scripted

## Scope

Reproduction starts from processed, analysis-ready data. Raw mass spectrometry files, the
Spectronaut processing that produced the normalized protein tables, and the FASTQ-to-count step of
the CRISPR screen are not included; they are available from the public repositories:

- mass spectrometry: MassIVE MSV000099551, ProteomeXchange PXD069731
- CRISPR screen sequencing: NCBI BioProject PRJNA1379893

## Licence

Code is MIT (`LICENSE`); data and generated outputs are CC BY 4.0 (`LICENSE-DATA`). `LICENSE-DATA`
states which paths each covers, and names the four external inputs that keep the terms of their
own sources.

## External data

Four inputs under `data/` come from outside this study. Each is reduced to what the analysis reads.

| file | source | version | note |
|---|---|---|---|
| `data/ibioid/go_string/*.tsv`, `data/ibioid/go_enrichment/*.tsv`, `data/ibiocofrac/m11/go_enrichment/*.tsv` | STRING, string-db.org | v12.0 | functional enrichment computed on this study's gene lists. STRING data are released under CC BY 4.0 |
| `data/external/syngo_1.3_Lrrk2_WT_ontologies_with_annotations.xlsx` | SynGO, syngoportal.org | release 1.3, run 2026-08-03 | synapse enrichment computed on this study's gene list |
| `data/external/PanglaoDB_markers_27_Mar_2020.tsv` | PanglaoDB, panglaodb.se | 27 March 2020 marker release | reduced to the mouse brain rows of the six cell types the marker construction uses, and to the four columns it reads |
| `data/external/annotation.BrainCellAtlas_Saunders_version_2018.04.01.csv` | Saunders et al. mouse brain atlas | 2018.04.01 | reduced to the substantia nigra rows |

`data/ibioid/orthogonal/recovered_partners_3DB.csv` records, per protein, whether BioGRID, IntAct
or STRING lists it as a physical interactor of the bait. Only those flags are included; no
interaction records or evidence text from those databases are reproduced here.

## Outputs

`CHECKSUMS.csv` records the size and md5 of every file the builders produce. SVG files carry a size
only: matplotlib writes a timestamp and randomised clip identifiers into SVG, so an SVG md5 is not
a stable check.
