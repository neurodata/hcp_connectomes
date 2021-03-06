from pathlib import Path

import boto3
from joblib import Parallel, delayed


def get_data(
    access_key_id,
    secret_access_key,
    output_path,
    prefix="HCP_1200/",
    n_jobs=1,
    verbose=False,
):
    """
    Do not hard code access key and secret key.

    Parameters
    ----------
    access_key_id : str

    secret_access_key : str

    output_path : str
        Path to the output files

    prefix : str
        One of {"HCP_1200/", "HCP_Retest/"}

    verbose : bool, default=False
    """
    bucket = "hcp-openaccess"

    def get_subjects():
        s3 = boto3.client(
            "s3",
            aws_access_key_id=access_key_id,
            aws_secret_access_key=secret_access_key,
        )
        continuation_token = None
        while True:
            list_kwargs = dict(Bucket=bucket, Prefix=prefix, Delimiter="/")
            if continuation_token:
                list_kwargs["ContinuationToken"] = continuation_token
            response = s3.list_objects_v2(**list_kwargs)
            yield from response.get("CommonPrefixes", [])
            if not response.get("IsTruncated"):  # At the end of the list?
                break
            continuation_token = response.get("NextContinuationToken")

    subject_list = [d["Prefix"] for d in get_subjects()]

    # Make output directories
    p = Path(output_path)

    def worker(sub_prefix):
        s3 = boto3.client(
            "s3",
            aws_access_key_id=access_key_id,
            aws_secret_access_key=secret_access_key,
        )
        sub_id = sub_prefix.split("/")[1]

        if verbose:
            print(f"Downloading Subject: {sub_id}...")

        # Get wmparc
        wmparc_key = s3.list_objects(
            Bucket="hcp-openaccess", Prefix=f"{sub_prefix}T1w/wmparc.nii.gz"
        )

        # Get eddy corrected files
        diffusion_keys = s3.list_objects(
            Bucket="hcp-openaccess", Prefix=f"{sub_prefix}T1w/Diffusion/", Delimiter="/"
        )

        try:
            to_download = [
                files["Key"]
                for files in wmparc_key["Contents"] + diffusion_keys["Contents"]
            ]
        except:
            return

        for key in to_download:
            filename = p / key.replace(prefix, "")
            if not filename.parents[0].exists():
                filename.parents[0].mkdir(parents=True)

            if verbose:
                print(f"Downloading File: {filename}...")

            s3.download_file(Bucket=bucket, Key=key, Filename=str(filename))

    if n_jobs == 1:
        for sub_prefix in subject_list:
            worker(sub_prefix)
    else:
        Parallel(n_jobs=n_jobs, verbose=1)(delayed(worker)(s) for s in subject_list)

