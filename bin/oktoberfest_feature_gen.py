#!/usr/bin/env python

"""
Generates features for PSM-rescoring using Oktoberfest.
The rescoring itself is suppresed by setting an unknown FDR estimation method.
"""

import argparse
import copy
import json
import logging
from pathlib import Path

import oktoberfest as ok
from oktoberfest import runner as ok_runner
import pandas as pd
import psm_utils
import psm_utils.io

OKTOBERFEST_UNKNOWN_FDR_ESTIMATION_METHOD = 'f{config.fdr_estimation_method} is not a valid rescoring tool, use either "percolator" or "mokapot"'
"""Oktoberfest's error message for unknown FDR estimation methods.
There is a typo in the oktberfest code, as it is not substituting the f-string correctly.
"""

def parse_str_bool(value: str) -> bool:
    """
    Parses a string argument to a boolean value.

    Arguments
    ---------
    value : str
        The string value to parse. Accepts 'true', 'false', '1', '

    Returns
    -------
    bool
        Returns True for 'true' or '1', and False for 'false' or '

    Raises
    ------
    ValueError
        If the value is not a valid boolean representation.
    """
    if value.lower() in ['true', '1']:
        return True
    elif value.lower() in ['false', '0']:
        return False
    else:
        raise ValueError(f"Invalid boolean value: {value}")

def argparse_setup() -> argparse.Namespace:
    """
    Creates the argument parser for the Oktoberfest feature generation script.
    """

    parser = argparse.ArgumentParser()
    # files
    parser.add_argument(
        "-psms-file", help="Input PSMs TSV file", required=True, type=Path
    )
    parser.add_argument(
        "-spectra-file",
        help="Corresponding spectrum file for PSMs file",
        required=True,
        type=Path,
    )

    # prediction
    parser.add_argument("-intensity-model", help="Koina intensity model", type=str)
    parser.add_argument("-irt-model", help="Koina IRT model", type=str)

    # mass spec parameters
    parser.add_argument(
        "-mass-tolerance",
        help="Defines the allowed tolerance between theoretical and experimentally observered fragment mass during peak annotation; default = 20 (FTMS), 40 (TOF), 0.35 (ITMS)",
        default=20,
        type=float,
    )
    parser.add_argument(
        "-mass-tolerance-unit",
        help="Defines the measure of tolerance, either “da” or “ppm”; default = da (mass analyzer is ITMS), ppm (mass analyzer is FTMS or TOF)",
        type=str,
        choices=["da", "ppm"],
    )

    parser.add_argument(
        "-is-timstof",
        help="If true, the spectra file type is set to 'd' (for timsTOF); otherwise, it defaults to 'mzml'",
        type=parse_str_bool
    )

    parser.add_argument(
        "-out-folder", help="Output folder for ", required=True, type=Path
    )

    return parser.parse_args()


def main():
    """
    Generates features for PSM-rescoring using Oktoberfest.
    The rescoring itself is suppresed by setting an unknown FDR estimation method.

    The resulting features can be found as `<output_folder>/none/rescore.tab`
    """

    args = argparse_setup()
    logging.basicConfig(level=logging.INFO)

    oktoberfest_input_csv_path = args.psms_file.with_suffix(".oktoberfest.input.csv")

    psms = psm_utils.io.read_file(args.psms_file)

    # Necessary columns according to the docs:
    # RAW_FILE,SCAN_NUMBER,MODIFIED_SEQUENCE,PRECURSOR_CHARGE,
    # SCAN_EVENT_NUMBER,MASS,SCORE,REVERSE,SEQUENCE,PEPTIDE_LENGTH
    oktoberfest_df = pd.DataFrame()

    # RAW_FILE,
    oktoberfest_df["RAW_FILE"] = [args.spectra_file.stem] * len(psms)

    # SCAN_NUMBER
    oktoberfest_df["SCAN_NUMBER"] = [psm.spectrum_id for psm in psms]

    # MODIFIED_SEQUENCE
    oktoberfest_df["MODIFIED_SEQUENCE"] = [
        psm.peptidoform.modified_sequence for psm in psms
    ]
    oktoberfest_df["MODIFIED_SEQUENCE"] = oktoberfest_df["MODIFIED_SEQUENCE"].apply(
        lambda x: x.replace("[UNIMOD:Carbamidomethyl]", "[UNIMOD:4]").replace(
            "[UNIMOD:Oxidation]", "[UNIMOD:35]"
        )
    )

    # PRECURSOR_CHARGE
    oktoberfest_df["PRECURSOR_CHARGE"] = [psm.get_precursor_charge() for psm in psms]

    # MASS
    oktoberfest_df["MASS"] = [psm.peptidoform.theoretical_mass for psm in psms]

    # SCORE
    # TODO: Does the psmutil score means higher is better?
    oktoberfest_df["SCORE"] = [psm.score for psm in psms]

    # REVERSE
    oktoberfest_df["REVERSE"] = [psm.is_decoy for psm in psms]

    # SEQUENCE
    oktoberfest_df["SEQUENCE"] = [psm.peptidoform.sequence for psm in psms]

    # PEPTIDE_LENGTH
    oktoberfest_df["PEPTIDE_LENGTH"] = [len(psm.peptidoform.sequence) for psm in psms]

    # free up some memory
    del psms

    psms_df = pd.read_csv(args.psms_file, sep="\t")

    # PROTEINS (not in the docs, but required by oktberfest)
    oktoberfest_df["PROTEINS"] = psms_df["protein_list"].apply(
        lambda x: x.replace("[", "").replace("]", "").replace("'", "")
    )  # remove the brackets and quotes

    # adding present rescoring columns
    present_rescoring_cols = list(
        filter(lambda x: x.startswith("rescoring:"), psms_df.columns)
    )
    for rescoring_col in present_rescoring_cols:
        oktoberfest_df[rescoring_col] = psms_df[rescoring_col]

    oktoberfest_df.to_csv(
        oktoberfest_input_csv_path,
        sep=",",
        index=False,
    )

    # free up some more memory as the dataframe is read fromt disk again
    del psms_df
    del oktoberfest_df

    # create the config file
    config_dict = copy.deepcopy(ok.utils.example_configs.RESCORING)

    # misc
    config_dict["output"] = str(args.out_folder)
    config_dict["num_threads"] = 1  # Set to 1 for debugging, can be increased later
    # mass spec parameters
    config_dict["mass_tolerance"] = args.mass_tolerance
    config_dict["unitMassTolerance"] = args.mass_tolerance_unit
    # predicition params
    config_dict["models"]["irt_model"] = args.irt_model
    config_dict["models"]["intensity_model"] = args.intensity_model
    # input params
    config_dict["inputs"]["search_results"] = str(oktoberfest_input_csv_path)
    config_dict["inputs"]["search_results_type"] = "Internal"
    config_dict["inputs"]["spectra"] = "./"
    config_dict["inputs"]["spectra_type"] = "d" if args.is_timstof else "mzml"
    # resocreing params
    # deliberately set to NONE, which will cause Oktoberfest to shut down before rescoring
    # by raising a ValueError which we can catch later.
    # This has the effect, that the generated features
    # are stored in the subfolder `results/none` of the output folder.
    config_dict["fdr_estimation_method"] = "NONE"
    config_dict["quantification"] = False
    config_dict["add_feature_cols"] = "all"

    config_path = Path("oktoberfest.config.json")

    with config_path.open("w", encoding="utf-8") as json_file:
        json_file.write(json.dumps(config_dict))

    try:
        ok_runner.run_job(config_path)
    except ValueError as e:
        # Catch the specific ValueError raised when "NONE" is used as fdr_estimation
        if str(e) != OKTOBERFEST_UNKNOWN_FDR_ESTIMATION_METHOD:
            raise e


if __name__ == "__main__":
    main()
