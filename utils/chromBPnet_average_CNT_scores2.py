#!/usr/bin/env python3
"""
average_contrib_scores.py

Average multiple h5 contribution-score files (counts and profile) into two averaged HDF5 outputs
(no_bias_model/{celltype}_averaged_contribution_scores_counts.h5 and
 no_bias_model/{celltype}_averaged_contribution_scores_profile.h5)

Designed to be memory-efficient (reads/writes in chunks), preserve attributes and group structure,
and produce compressed output (gzip+shuffle) to reduce file size.

Usage:
    python average_contrib_scores.py /full/path/to/chromBPct_dir CELLTYPE [--dtype float32] [--compress-level 6]

Author: ChatGPT (reworked for Bruno's workflow)
"""
import os
import sys
import argparse
import math
import h5py
import numpy as np
from pathlib import Path
import logging

# --- Config / helpers ---
DEFAULT_COMPRESS_LEVEL = 6
# target bytes to read per chunk (tuneable). ~100MB is a reasonable default.
TARGET_BYTES_PER_CHUNK = 100 * 1024**2

logging.basicConfig(format="%(asctime)s %(levelname)s: %(message)s", level=logging.INFO)


def parse_args():
    p = argparse.ArgumentParser(description="Memory-friendly average of HDF5 contribution score files.")
    p.add_argument("chromBPct_dir", help="Top-level directory (your chromBPct dir).")
    p.add_argument("celltype", help="Celltype (used to build output filenames).")
    p.add_argument("--dtype", choices=("float32", "float64"), default="float32",
                   help="Output dtype for averaged datasets (default: float32).")
    p.add_argument("--compress-level", type=int, default=DEFAULT_COMPRESS_LEVEL,
                   help="Gzip compression level for output (0-9).")
    return p.parse_args()


def list_expected_files(chromBPct_dir, prefix="no_bias_model"):
    base = Path(chromBPct_dir) / prefix
    counts_files = [base / f"fold_{i}" / "contribution_scores.counts_scores.h5" for i in range(5)]
    profile_files = [base / f"fold_{i}" / "contribution_scores.profile_scores.h5" for i in range(5)]
    return counts_files, profile_files


def check_files_exist(files):
    paths = [Path(f) for f in files]
    missing = [str(p) for p in paths if not p.exists()]
    if missing:
        logging.error("Missing files:\n" + "\n".join(missing))
        return False
    return True

def choose_chunk_size(shape, dtype, target_bytes=TARGET_BYTES_PER_CHUNK):
    """
    Choose a chunk along axis 0 so that approx target_bytes are read per block.
    Returns an integer chunk size for axis0 (>=1).
    """
    if len(shape) == 0:
        return 1
    itemsize = np.dtype(dtype).itemsize
    # If dataset has axis-0 length small, return full length
    axis0 = shape[0]
    # bytes per one slice of axis0:
    bytes_per_slice = itemsize
    for s in shape[1:]:
        bytes_per_slice *= int(s)
    if bytes_per_slice == 0:
        return axis0 or 1
    chunk = max(1, int(target_bytes // bytes_per_slice))
    return min(axis0, chunk)


def average_datasets(input_paths, group_key, dset_key, out_dset, out_dtype):
    """
    Compute mean across input_paths for dataset located at group_key/dset_key.
    Writes into out_dset (an h5py.Dataset already created).
    """
    # open first to get shape/dtype safely
    with h5py.File(input_paths[0], "r") as fh0:
        ds0 = fh0[group_key][dset_key]
        shape = ds0.shape
        src_dtype = ds0.dtype

    # handle scalars (0-d)
    if len(shape) == 0:
        # read scalar from each and compute mean
        vals = []
        for p in input_paths:
            with h5py.File(p, "r") as fh:
                vals.append(np.array(fh[group_key][dset_key][()]))
        meanval = np.mean(np.array(vals), axis=0).astype(out_dtype)
        out_dset[()] = meanval
        return

    # choose chunk along axis 0
    chunk0 = choose_chunk_size(shape, src_dtype, TARGET_BYTES_PER_CHUNK)
    logging.debug(f"dataset {group_key}/{dset_key} shape={shape} chunk0={chunk0}")

    # accumulator in float64 for numeric stability
    accum_dtype = np.float64
    accum = np.zeros((min(chunk0, shape[0]),) + tuple(shape[1:]), dtype=accum_dtype)  # reused buffer

    # iterate slices along axis0
    total0 = shape[0]
    start = 0
    while start < total0:
        this_chunk = min(chunk0, total0 - start)
        accum[:this_chunk] = 0.0
        for p in input_paths:
            with h5py.File(p, "r") as fh:
                src = fh[group_key][dset_key]
                sl = tuple([slice(start, start + this_chunk)] + [slice(None)] * (len(shape) - 1))
                # read to float64
                accum_slice = np.array(src[sl], dtype=accum_dtype)
                accum[:this_chunk] += accum_slice
        # mean
        mean_block = (accum[:this_chunk] / float(len(input_paths))).astype(out_dtype)
        # write to out_dset
        dest_sl = tuple([slice(start, start + this_chunk)] + [slice(None)] * (len(shape) - 1))
        out_dset[dest_sl] = mean_block
        start += this_chunk


def process_file_list(input_paths, output_path, out_dtype_str="float32", compress_level=6):
    input_paths = [Path(p) for p in input_paths]
    out_dtype = np.dtype(out_dtype_str)
    out_path = Path(output_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    # quick checks
    if not check_files_exist(input_paths):
        raise FileNotFoundError("One or more input files are missing. Aborting.")

    # open reference to iterate groups/datasets
    with h5py.File(input_paths[0], "r") as ref:
        groups = list(ref.keys())

    logging.info(f"Processing {len(input_paths)} files -> {out_path}. Groups found: {groups}")

    # Create output file
    with h5py.File(out_path, "w") as out_f:
        # copy top-level attrs if any
        with h5py.File(input_paths[0], "r") as ref:
            for k, v in ref.attrs.items():
                out_f.attrs[k] = v

        for group_key in groups:
            logging.info(f"Working on group: {group_key}")
            # ensure group exists in output
            out_group = out_f.create_group(group_key)
            # copy group attrs
            with h5py.File(input_paths[0], "r") as ref:
                for k, v in ref[group_key].attrs.items():
                    out_group.attrs[k] = v

            # iterate dataset keys in this group (use reference file)
            with h5py.File(input_paths[0], "r") as ref:
                ds_keys = list(ref[group_key].keys())

            for dset_key in ds_keys:
                logging.info(f"  Averaging dataset: {dset_key}")
                # get shape and create dataset in output with compression
                with h5py.File(input_paths[0], "r") as ref:
                    src = ref[group_key][dset_key]
                    shape = src.shape
                    attrs = dict(src.attrs)
                # choose chunks for output (let h5py pick automatic chunk dims by None)
                # create dataset
                create_kwargs = dict(dtype=out_dtype, compression="gzip",
                                     compression_opts=compress_level, shuffle=True)
                out_dset = out_group.create_dataset(dset_key, shape=shape, **create_kwargs)
                # copy dataset attrs (if any) after creation
                for k, v in attrs.items():
                    try:
                        out_dset.attrs[k] = v
                    except Exception:
                        logging.debug(f"Could not copy attribute {k} for {group_key}/{dset_key}")

                # compute average (chunked)
                try:
                    average_datasets(input_paths, group_key, dset_key, out_dset, out_dtype)
                except Exception as e:
                    logging.error(f"Error averaging dataset {group_key}/{dset_key}: {e}")
                    raise

    logging.info(f"Finished writing averaged file: {out_path}")


def main():
    args = parse_args()
    chromBPct_dir = args.chromBPct_dir
    celltype = args.celltype
    out_dtype = args.dtype
    compress_level = args.compress_level

    counts_files, profile_files = list_expected_files(chromBPct_dir)
    counts_out = Path(chromBPct_dir) / "no_bias_model" / f"{celltype}_averaged_contribution_scores_counts.h5"
    profile_out = Path(chromBPct_dir) / "no_bias_model" / f"{celltype}_averaged_contribution_scores_profile.h5"

    logging.info("Counts files:\n  " + "\n  ".join(str(p) for p in counts_files))
    logging.info("Profile files:\n  " + "\n  ".join(str(p) for p in profile_files))

    # process counts
    process_file_list(counts_files, counts_out, out_dtype_str=out_dtype, compress_level=compress_level)
    # process profile
    process_file_list(profile_files, profile_out, out_dtype_str=out_dtype, compress_level=compress_level)


if __name__ == "__main__":
    main()
