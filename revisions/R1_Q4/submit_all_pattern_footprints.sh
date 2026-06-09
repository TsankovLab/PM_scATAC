#!/bin/bash
# submit_all_pattern_footprints.sh  (R1_Q4)
# Submits one GPU job per KLRC1_NK counts-modisco pattern.
# Skips patterns whose output directory already has all 30 fold H5s.

SCRIPT=/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/git_repo_claude/R1_Q4/run_footprint_generic.sh
MESO_DIR=/sc/arion/scratch/giottb01/Xmen/meso
LOG_DIR=/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/git_repo_claude/logs

chmod +x "${SCRIPT}"
mkdir -p "${LOG_DIR}"

# pattern_id  core_fwd        core_rc
PATTERNS=(
"pos_patterns.pattern_0  ACTTCCTG CAGGAAGT"
"pos_patterns.pattern_1  AGGGGGCG CGCCCCCT"
"pos_patterns.pattern_2  AAACCACA TGTGGTTT"
"pos_patterns.pattern_3  GGGCGGGG CCCCGCCC"
"pos_patterns.pattern_4  TTTCACTT AAGTGAAA"
"pos_patterns.pattern_5  TGAGTCAT ATGACTCA"
"pos_patterns.pattern_6  TGACGTAA TTACGTCA"
"pos_patterns.pattern_7  TGATTGGC GCCAATCA"
"pos_patterns.pattern_8  CGCAGGCG CGCCTGCG"
"pos_patterns.pattern_9  AGGTGTGA TCACACCT"
"pos_patterns.pattern_10 GGAAATCC GGATTTCC"
"pos_patterns.pattern_11 CCACATCC GGATGTGG"
"pos_patterns.pattern_12 GGGAATTG CAATTCCC"
"pos_patterns.pattern_13 AGATGGCG CGCCATCT"
"pos_patterns.pattern_14 TTTGGTTT AAACCAAA"
"pos_patterns.pattern_15 GAAACCAC GTGGTTTC"
"pos_patterns.pattern_16 ATGGCAAC GTTGCCAT"
"pos_patterns.pattern_17 GTCACGTG CACGTGAC"
"pos_patterns.pattern_18 TTCCTTCC GGAAGGAA"
"pos_patterns.pattern_19 CCACAAAC GTTTGTGG"
"pos_patterns.pattern_20 TGACATCA TGATGTCA"
"pos_patterns.pattern_21 TCCTGTGG CCACAGGA"
"pos_patterns.pattern_22 CAAACCAC GTGGTTTG"
"pos_patterns.pattern_23 TCTCGCGA TCGCGAGA"
"pos_patterns.pattern_24 TGAGTCAT ATGACTCA"
"pos_patterns.pattern_25 TTTTTTTT AAAAAAAA"
"pos_patterns.pattern_26 TAACCACA TGTGGTTA"
"pos_patterns.pattern_27 GGGTGTGA TCACACCC"
"pos_patterns.pattern_28 CGCTTCCG CGGAAGCG"
"pos_patterns.pattern_29 TTCTGGGA TCCCAGAA"
"pos_patterns.pattern_30 CCACATCC GGATGTGG"
"pos_patterns.pattern_31 AAACCACA TGTGGTTT"
"pos_patterns.pattern_32 AGGTGTGA TCACACCT"
"pos_patterns.pattern_33 TGTGGTTT AAACCACA"
"pos_patterns.pattern_34 GGAAACCC GGGTTTCC"
"pos_patterns.pattern_35 ATTTCCTT AAGGAAAT"
"pos_patterns.pattern_36 TGTCACTT AAGTGACA"
"pos_patterns.pattern_37 TGAGTCAT ATGACTCA"
"pos_patterns.pattern_38 TCCCCGCC GGCGGGGA"
"pos_patterns.pattern_39 GGTGTGAA TTCACACC"
"pos_patterns.pattern_40 GGTGTGAC GTCACACC"
"pos_patterns.pattern_41 GGATGTGT ACACATCC"
"neg_patterns.pattern_0  CTTGGAAT ATTCCAAG"
"neg_patterns.pattern_1  ATTCCATG CATGGAAT"
"neg_patterns.pattern_2  CTGGGATT AATCCCAG"
"neg_patterns.pattern_3  CACATCCT AGGATGTG"
)

for entry in "${PATTERNS[@]}"; do
    read -r pat_id core_fwd core_rc <<< "${entry}"
    pat_safe=$(echo "${pat_id}" | tr './' '__')
    out_dir="${MESO_DIR}/footprints_${pat_safe}"

    # Skip if all 30 fold H5s already exist
    n_done=$(ls "${out_dir}"/*_footprints.h5 2>/dev/null | grep -v average | wc -l)
    if [ "${n_done}" -ge 30 ]; then
        echo "SKIP ${pat_id} (${n_done}/30 folds done)"
        continue
    fi

    job_name="fp_${pat_safe:0:18}"
    echo "SUBMIT ${pat_id}  core=${core_fwd}  done=${n_done}/30"

    bsub \
        -J  "${job_name}" \
        -P  acc_Tsankov_Normal_Lung \
        -q  gpu \
        -n  4 \
        -W  12:00 \
        -gpu "num=1" \
        -R  "a100" \
        -R  "rusage[mem=32000]" \
        -R  "span[hosts=1]" \
        -o  "${LOG_DIR}/fp_${pat_safe}_%J.out" \
        -e  "${LOG_DIR}/fp_${pat_safe}_%J.err" \
        "${SCRIPT}" "${pat_id}" "${core_fwd}" "${core_rc}"
done

echo "=== Submission complete ==="
