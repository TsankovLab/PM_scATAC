# R2_Q3 — subclonal structure in the scATAC tumours (epiAneufinder)

Reviewer 2, Question 3. epiAneufinder subclones in every scATAC tumour, TF motifs that differ between the subclones, and
spatial validation of the P4 chr8q clone in Visium. Method text: `METHODS_subclones.md`.

## How to reproduce

```bash
cd .../git_repo_claude/R2_Q3
./run_all.sh          # steps 3, 4, 9 (three variants), 10 (two variants), 11
```

Run on a node with >= 8 GB (step 3 holds an 11,331 x 11,331 distance matrix for P23). `run_all.sh` starts at step 3
because `out_5Mb/` holds the epiAneufinder output; to redo steps 1-2:

```bash
./run_step.sh 01_extract_fragments.R          # ~10 min, writes ~11 GB into frags/
for S in P1 P4 P5 P8 P10 P11 P12 P14 P23; do
  bsub -P acc_Tsankov_Normal_Lung -q premium -n 8 -R "rusage[mem=8000] span[hosts=1]" \
       -W 6:00 -J epi_$S -o logs/epi_$S.out -e logs/epi_$S.err \
       "$PWD/run_step.sh 02_run_epianeufinder.R $S"
done
```

## Scripts

| step | script | what it does | output |
|---|---|---|---|
| — | `00_common.R` | paths, epiAneufinder parameters, clone-split function, palette; sourced by every step | — |
| 1 | `01_extract_fragments.R` | per-cell fragments of each tumour from `tumor_compartment/scatac_ArchR` | `frags/` |
| 2 | `02_run_epianeufinder.R <S>` | epiAneufinder at 5 Mb, one sample | `out_5Mb/<S>/` |
| 3 | `03_clone_calls.R` | defines the clones (primary split accepted when an arm separates by >= 0.55) | `epi_clone_labels.csv`, `epi_clone_decision.csv`, `epi_clone_profiles.rds`, ... |
| 4 | `04_figure_subclones_overview.R` | CNV landscape of all clones and what separates each pair | `Plots/epiAneufinder_subclones_overview.pdf` |
| 9 | `09_P4_visium_chr8q.R` | P4 chr8q clone on the four Visium sections; `MAL_CALL=label` uses the label-transfer annotation, `SCORE=relative` the genome-relative chr8q score | `Plots/R2_Q3_P4_visium_chr8q{,_label,_label_relative}.pdf` |
| 10 | `10_P4_visium_chr8q_clustering.R` | are chr8q-amp spots neighbours among malignant spots? (999 label shuffles) | `Plots/R2_Q3_P4_visium_chr8q_clustering{,_label}.pdf/.png` |
| 11 | `11_tf_subclone_wilcoxon.R` | per tumour, Wilcoxon of chromVAR z between the two clones; motifs ranked by the number of tumours with FDR < 0.05; TF must be expressed in the tumour's malignant scRNA | `Plots/epi_TF_variability_wilcoxon.pdf/.png` |

`run_step.sh` activates the `meso_scatac` environment and runs one script.
