from pathlib import Path

import boto3
fromb joblib import Parallel, delayed

def get_data(access_key_id, secret_access_key, output_path, n_jobs=1, verbose=False):
    """
    Do not hard code access key and secret key.

    Parameters
    ----------
    access_key_id : str

    secret_access_key : str

    output_path : str
        Path to the output files

    verbose : bool, default=False
    """
    bucket = "hcp-openaccess"
    prefix = "HCP_1200/"

    s3 = boto3.client(
        "s3", aws_access_key_id=access_key_id, aws_secret_access_key=secret_access_key
    )
    result = s3.list_objects(Bucket=bucket, Prefix=prefix, Delimiter="/")

    subject_list = []
    for o in result.get("CommonPrefixes"):
        subject_list.append(o.get("Prefix"))

    # Make output directories
    p = Path(output_path)

    def worker(sub_prefix):
        sub_id = sub_prefix.split("/")[1]

        if verbose:
            print(f"Downloading Subject: {sub_id}...")

        sub_path = p / sub_id
        sub_path.mkdir(parents=True, exist_ok=True)

        # Get wmparc
        wmparc_key = s3.list_objects(
            Bucket="hcp-openaccess", Prefix=f"{sub_prefix}T1w/wmparc.nii.gz"
        )

        # Get eddy corrected files
        diffusion_keys = s3.list_objects(
            Bucket="hcp-openaccess", Prefix=f"{sub_prefix}T1w/Diffusion/", Delimiter="/"
        )

        to_download = [
            files["Key"]
            for files in wmparc_key["Contents"] + diffusion_keys["Contents"]
        ]

        for key in to_download:
            filename = p / key.replace(prefix, "")
            if not filename.parents[0].exists():
                filename.parents[0].mkdir()

            if verbose:
                print(f"Downloading File: {filename}...")

            s3.download_file(Bucket=bucket, Key=key, Filename=filename)

    if n_jobs == 1:
        for sub_prefix in subject_list:
            worker(sub_prefix)

    else:
        Parallel(n_jobs=n_jobs)(delayed(worker)(sub_prefix) for sub_prefix in subject_list)