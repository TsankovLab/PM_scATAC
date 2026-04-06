#!/usr/bin/env python3
"""
average_contrib_scores_tables.py

Memory-efficient averaging of multiple ChromBPNet HDF5 contribution-score files
(counts and profile) using PyTables (Blosc-compatible).

Usage:
    python average_contrib_scores_tables.py /full/path/to/chromBPct_dir CELLTYPE [--dtype float32] [--chunk-bytes 100000000]

Author: ChatGPT reworked for Bruno's workflow
"""
import os
import argparse
import logging
import tables as tb
import numpy as np
from pathlib import Path

logging.basicConfig(format="%(asctime)s %(levelname)s: %(message)s", level=logging.INFO)

# --- helpers ---
DEFAULT_CHUNK_BYTES = 100 * 1024**2  # 100 MB


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("chromBPct_dir", help="Top-level chromBPnet directory")
    p.add_argument("celltype", help="Celltype (used to build output filenames)")
    p.add_argument("--dtype", choices=("float32", "float64"), default="float32")
    p.add_argument("--chunk-bytes", type=int, default=DEFAULT_CHUNK_BYTES)
    return p.parse_args()


def list_expected_files(chromBPct_dir):
    base = Path(chromBPct_dir) / "no_bias_model"
    counts_files = [base / f"fold_{i}" / "contribution_scores.counts_scores.h5" for i in range(5)]
    profile_files = [base / f"fold_{i}" / "contribution_scores.profile_scores.h5" for i in range(5)]
    return counts_files, profile_files


def check_files_exist(files):
    missing = [str(f) for f in files if not Path(f).exists()]
    if missing:
        logging.error("Missing files:\n" + "\n".join(missing))
        return False
    return True


def choose_chunk_size(shape, dtype, target_bytes=DEFAULT_CHUNK_BYTES):
    """Return axis0 chunk size to approximate target_bytes per read."""
    if len(shape) == 0:
        return 1
    itemsize = np.dtype(dtype).itemsize
    bytes_per_slice = itemsize
    for s in shape[1:]:
        bytes_per_slice *= int(s)
    if bytes_per_slice == 0:
        return shape[0] or 1
    return max(1, min(shape[0], target_bytes // bytes_per_slice))


def average_datasets(input_paths, group_name, dset_name, out_node, out_dtype, chunk_bytes):
    """Average dataset across multiple input files in chunks."""
    # get shape from first file
    with tb.open_file(input_paths[0], "r") as f:
        ds0 = f.get_node(f"/{group_name}/{dset_name}")
        shape = ds0.shape

    chunk0 = choose_chunk_size(shape, out_dtype, chunk_bytes)
    accum_dtype = np.float64
    accum = np.zeros((min(chunk0, shape[0]),) + tuple(shape[1:]), dtype=accum_dtype)

    total0 = shape[0]
    start = 0
    while start < total0:
        this_chunk = min(chunk0, total0 - start)
        accum[:this_chunk] = 0.0
        for p in input_paths:
            with tb.open_file(p, "r") as f:
                src = f.get_node(f"/{group_name}/{dset_name}")
                sl = slice(start, start + this_chunk)
                data = src[sl]
                accum[:this_chunk] += data.astype(accum_dtype)
        out_node[start:start+this_chunk] = (accum[:this_chunk] / len(input_paths)).astype(out_dtype)
        start += this_chunk


def process_file_list(input_paths, output_path, dtype_str, chunk_bytes):
    out_dtype = np.dtype(dtype_str)
    input_paths = [str(p) for p in input_paths]
    if not check_files_exist(input_paths):
        raise FileNotFoundError("Missing input files. Aborting.")

    # get group names (store strings only)
    with tb.open_file(input_paths[0], "r") as ref:
        groups = [g._v_name for g in ref.root._f_list_nodes()]

    logging.info(f"Processing {len(input_paths)} files -> {output_path}. Groups: {groups}")

    # open output file
    with tb.open_file(output_path, "w") as out_f:
        for group_name in groups:
            logging.info(f"Working on group: {group_name}")
            out_group = out_f.create_group("/", group_name)
            # copy group attrs
            with tb.open_file(input_paths[0], "r") as ref:
                group_node = ref.get_node(f"/{group_name}")
                for k in group_node._v_attrs._f_list():
                    out_group._v_attrs[k] = group_node._v_attrs[k]

                # iterate datasets
                for dset_node in group_node._f_list_nodes():
                    dset_name = dset_node._v_name
                    logging.info(f"  Averaging dataset: {dset_name}")
                    # create output array
                    arr_shape = dset_node.shape
                    out_arr = out_f.create_carray(out_group, dset_name,
                                                  atom=tb.Atom.from_dtype(out_dtype),
                                                  shape=arr_shape,
                                                  filters=tb.Filters(complevel=6, complib="blosc"),
                                                  chunkshape=None)
                    average_datasets(input_paths, group_name, dset_name, out_arr, out_dtype, chunk_bytes)

    logging.info(f"Finished writing averaged file: {output_path}")


def main():
    args = parse_args()
    chromBPct_dir = args.chromBPct_dir
    celltype = args.celltype
    dtype = args.dtype
    chunk_bytes = args.chunk_bytes

    counts_files, profile_files = list_expected_files(chromBPct_dir)

    counts_out = Path(chromBPct_dir) / "no_bias_model" / f"{celltype}_averaged_contribution_scores_counts.h5"
    profile_out = Path(chromBPct_dir) / "no_bias_model" / f"{celltype}_averaged_contribution_scores_profile.h5"

    logging.info("Counts files:\n  " + "\n  ".join(str(p) for p in counts_files))
    logging.info("Profile files:\n  " + "\n  ".join(str(p) for p in profile_files))

    process_file_list(counts_files, counts_out, dtype, chunk_bytes)
    process_file_list(profile_files, profile_out, dtype, chunk_bytes)


if __name__ == "__main__":
    main()
