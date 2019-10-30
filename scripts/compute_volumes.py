from argparse import ArgumentParser

from pathlib import Path
import os

import nibabel as nib
import numpy as np
import pandas as pd

from joblib import Parallel, delayed


def compute_brain_volumes(
    input_path, output_path, parcellation_file, n_jobs=-2, verbose=1
):
    """
    Tissue mask values
    white = 3
    gray = 2
    csf = 1
    """
    input_path = Path(input_path)
    output_path = Path(output_path)
    parcellation_path = Path(parcellation_file)

    if not output_path.is_dir():
        output_path.mkdir(parents=True)

    parcellation_img = nib.load(str(parcellation_path)).get_fdata()

    subjects = [x.name.split("-")[-1] for x in sorted(list(input_path.glob("*sub*")))]

    def compute_per_subject(subject):
        mask_path = input_path / f"sub-{subject}/masks/sub-{subject}_tissue_mask.nii.gz"
        tissue_mask = nib.load(str(mask_path)).get_fdata()

        sizes = []
        for roi in np.unique(parcellation_img)[1:]:
            tmp = (parcellation_img == roi) * (tissue_mask >= 2)
            sizes.append(tmp.sum())

        return sizes

    res = Parallel(n_jobs=n_jobs, verbose=verbose)(
        delayed(compute_per_subject)(subject) for subject in subjects
    )

    df = pd.DataFrame(
        np.array(res),
        index=subjects,
        columns=np.unique(parcellation_img)[1:].astype(int),
    )

    output_file = output_path / f"{parcellation_path.name.split('.nii')[0]}.csv"
    df.to_csv(output_file)


def main(input_path, output_path, parcellation_path):
    parcellations = Path(parcellation_path)

    all_parcellations = sorted(list(parcellations.glob("*1x1x1.nii.gz*")))
    exclude_list = ["DS", "yeo", "tissue", "slab", "hemispheric"]
    parcellation_names = [
        str(x.name)
        for x in all_parcellations
        if not any(substring in str(x) for substring in exclude_list)
    ]

    for parcel in parcellation_names:
        print(f"\nComputing brain volumes for {parcel.split('.')[0]}")
        parcellation_name = parcellations / parcel
        compute_brain_volumes(input_path, output_path, parcellation_name)


if __name__ == "__main__":
    parser = ArgumentParser(description="This is a script for computing brain volumes.")
    parser.add_argument(
        "input_dir",
        help="The directory with the input dataset"
        " formatted according to the BIDS standard.",
    )
    parser.add_argument(
        "output_dir",
        help="The directory where the output "
        "files should be stored. If you are running group "
        "level analysis this folder should be prepopulated "
        "with the results of the participant level analysis.",
    )
    parser.add_argument(
        "parcellation_dir", help="The neuroparc directoary containing parcellations."
    )

    result = parser.parse_args()
    inDir = result.input_dir
    outDir = result.output_dir
    parcDir = result.parcellation_dir

    main(inDir, outDir, parcDir)
