#!/usr/bin/env python

# imports
import argparse
import pandas as pd
from psm_utils.io import read_file

def argparse_setup():
    parser = argparse.ArgumentParser()
    parser.add_argument("-in_file", help="Input file name")
    parser.add_argument("-out_file", help="Output file name")
    parser.add_argument("-searchengine", help="The searchengine (only mandatory for comet, maxquant, msgfplus)")

    return parser.parse_args()


if __name__ == "__main__":
    args = argparse_setup()

    input_file = args.in_file
    output_file = args.out_file
    searchengine = args.searchengine

    # read in psm_utils TSV file
    psm_list = read_file(input_file, filetype="tsv")

    # calculate q-values using score (higher score better)
    psm_list.calculate_qvalues(reverse=True)
    psm_list.set_ranks()

    # filter by q-value, rank and target-decoy
    psm_df = psm_list.to_dataframe()
    psm_df = psm_df[(psm_df["qvalue"] < 0.01) & (psm_df["rank"] <= 1) & (psm_df["is_decoy"] == False)]

    # we need only these few data in the following steps
    small_df = pd.DataFrame(psm_df["peptidoform"].astype(str))
    small_df["spectrum_id"] = psm_df["spectrum_id"]
    small_df["qvalue"] = psm_df["qvalue"]

    # change the spectrum_id to the pure number (some searchengines use the format directly from mzML)
    if searchengine in {"msamanda", "msgfplus", "sage", "xtandem"}:
        small_df["spectrum_id"] = small_df["spectrum_id"].str.extract(r'.*scan=(\d+)')


    # write out the dataframe to file
    small_df.to_csv(output_file, sep="\t", index=False)
