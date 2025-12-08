#!/usr/bin/env python3
"""
average_contrib_scores_tables.py

Average multiple HDF5 contribution-score files (counts and profile) into two averaged HDF5 outputs
using PyTables (Blosc-compressed).

Usage:
    python average_contrib_scores_tables.py /full/path/to/chromBPct_dir CELLTYPE [--dtype float32] [--chunk-bytes 100MB]

Author: ChatGPT (adapted for Bruno's workflow)
"""
import os
import sys
import argparse
import logging
import numpy as np
import tables as tb
from pathlib import Path

logging.basicConfig(format="%(asctime)s %(levelname)s: %(message)s", level=logging.INFO)

# --- defaults ---
DEFAULT_CHUNK_BYTES = 100 * 1024**2  # ~100 MB

def parse_args():
    p = argparse.ArgumentParser(description="Average HDF5 contribution-score files using PyTables (Blosc).")
    p.add_argument("chromBPct_dir", help="Top-level directory for chromBPct.")
    p.add_argument("celltype", help="Celltype for output filenames.")
    p.add_argument("--dtype", choices=("float32", "float64"), default="float32", help="Output dtype.")
    p.add_argument("--chunk-bytes", type=int, default=DEFAULT_CHUNK_BYTES, help="Approx bytes per chunk for averaging.")
    return p.parse_args()

def list_expected_files(chromBPct_dir):
    base = Path(chromBPct_dir) / "no_bias_model"
    counts_files = [base / f"fold_{i}" / "contribution_scores.counts_scores.h5" for i in range(5)]
    profile_files = [base / f"fold_{i}" / "contribution_scores.profile_scores.h5" for i in range(5)]
    return counts_files, profile_files

def check_files_exist(files):
    missing = [str(f) for f in files if not f.exists()]
    if missing:
        logging.error("Missing files:\n" + "\n".join(missing))
        return False
    return True

def choose_chunk_size(shape, dtype, target_bytes):
    itemsize = np.dtype(dtype).itemsize
    bytes_per_slice = itemsize
    for s in shape[1:]:
        bytes_per_slice *= s
    if bytes_per_slice == 0:
        return shape[0]
    chunk = max(1, int(target_bytes // bytes_per_slice))
    return min(shape[0], chunk)

def average_datasets(input_paths, group, dset_key, out_node, out_dtype, chunk_bytes):
    # Open first file to get shape, compression info
    with tb.open_file(input_paths[0], "r") as f:
        node = f.get_node(f"/{group}/{dset_key}")
        shape = node.shape

    chunk0 = choose_chunk_size(shape, node.dtype, chunk_bytes)
    total0 = shape[0]

    logging.info(f"Dataset {group}/{dset_key} shape={shape}, chunk0={chunk0}")

    start = 0
    while start < total0:
        this_chunk = min(chunk0, total0 - start)
        accum = np.zeros((this_chunk, *shape[1:]), dtype=np.float64)
        for path in input_paths:
            with tb.open_file(path, "r") as f_in:
                arr = f_in.get_node(f"/{group}/{dset_key}")[start:start+this_chunk]
                accum[:arr.shape[0]] += arr
        out_node[start:start+this_chunk] = (accum / len(input_paths)).astype(out_dtype)
        start += this_chunk

def process_file_list(input_paths, output_path, out_dtype="float32", chunk_bytes=DEFAULT_CHUNK_BYTES):
    input_paths = [Path(p) for p in input_paths]
    if not check_files_exist(input_paths):
        raise FileNotFoundError("Missing input files.")

    out_dtype_np = np.dtype(out_dtype)
    out_path = Path(output_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    # Open first file to get groups/datasets
    with tb.open_file(input_paths[0], "r") as ref:
        groups = list(ref.root._f_list_nodes())

    logging.info(f"Processing {len(input_paths)} files -> {out_path}. Groups: {[g._v_name for g in groups]}")

    with tb.open_file(out_path, "w") as out_f:
        for group_node in groups:
            group_name = group_node._v_name
            logging.info(f"Working on group: {group_name}")
            out_group = out_f.create_group("/", group_name)
            # copy group attrs
            for k in group_node._v_attrs._f_list():
                out_group._v_attrs[k] = group_node._v_attrs[k]

            # iterate datasets
            for dset_node in group_node._f_list_nodes():
                dset_name = dset_node._v_name
                logging.info(f"  Averaging dataset: {dset_name}")

                # preserve compression and filters
                fobj = dset_node.filters
                filters = tb.Filters(complevel=fobj.complevel if fobj else 0,
                                     complib=fobj.complib if fobj else "blosc",
                                     shuffle=fobj.shuffle if fobj else True)

                atom = tb.Atom.from_dtype(out_dtype_np)
                out_node = out_f.create_carray(out_group, dset_name, atom=atom,
                                               shape=dset_node.shape, filters=filters)
                # copy attrs
                for k in dset_node._v_attrs._f_list():
                    out_node._v_attrs[k] = dset_node._v_attrs[k]

                # average in chunks
                average_datasets(input_paths, group_name, dset_name, out_node, out_dtype, chunk_bytes)

    logging.info(f"Finished writing: {out_path}")

def main():
    args = parse_args()
    chromBPct_dir = args.chromBPct_dir
    celltype = args.celltype
    out_dtype = args.dtype
    chunk_bytes = args.chunk_bytes

    counts_files, profile_files = list_expected_files(chromBPct_dir)
    counts_out = Path(chromBPct_dir) / "no_bias_model" / f"{celltype}_averaged_contribution_scores_counts.h5"
    profile_out = Path(chromBPct_dir) / "no_bias_model" / f"{celltype}_averaged_contribution_scores_profile.h5"

    logging.info("Counts files:\n" + "\n".join(str(f) for f in counts_files))
    logging.info("Profile files:\n" + "\n".join(str(f) for f in profile_files))

    process_file_list(counts_files, counts_out, out_dtype, chunk_bytes)
    process_file_list(profile_files, profile_out, out_dtype, chunk_bytes)

if __name__ == "__main__":
    main()
