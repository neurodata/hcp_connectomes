import os
from pathlib import Path

import nibabel as nib
from ndmg.utils import reg_utils as mgru


def match_target_vox_res(img_file, vox_size="1mm", sens="t1w"):
    """Reslices input MRI file if it does not match the targeted voxel resolution. Can take dwi or t1w scans.
    
    Parameters
    ----------
    img_file : str
        path to file to be resliced
    vox_size : str
        target voxel resolution ('2mm' or '1mm')
    namer : name_resource
        name_resource variable containing relevant directory tree information
    sens : str
        type of data being analyzed ('dwi' or 'func')
    
    Returns
    -------
    str
        location of potentially resliced image
    """
    from dipy.align.reslice import reslice

    # Check dimensions
    img = nib.load(img_file)
    data = img.get_fdata()
    affine = img.affine
    hdr = img.header
    zooms = hdr.get_zooms()[:3]
    if vox_size == "1mm":
        new_zooms = (1.0, 1.0, 1.0)
    elif vox_size == "2mm":
        new_zooms = (2.0, 2.0, 2.0)

    if (abs(zooms[0]), abs(zooms[1]), abs(zooms[2])) != new_zooms:
        # print("Reslicing image " + img_file + " to " + vox_size + "...")

        data2, affine2 = reslice(data, affine, zooms, new_zooms)
        img2 = nib.Nifti1Image(data2, affine=affine2)
        nib.save(img2, img_file)
    else:
        nib.save(img, img_file)


def register_t1w_2_mni(
    input_path, output_path, subject, ses=1, nonlinear=False, vox_size="1mm"
):
    """
    Parameters
    ----------
    input_path : str
        directory to bids
    output_path : str
        output bids
    """
    # deal with paths
    input_path = Path(input_path)
    output_path = Path(output_path) / f"sub-{subject}"

    if not output_path.is_dir():
        output_path.mkdir(parents=True)

    FSLDIR = os.environ["FSLDIR"]

    # files names
    input_t1w = (
        input_path
        / f"sub-{subject}"
        / f"ses-{ses}"
        / "anat"
        / f"sub-{subject}_ses-{ses}_T1w.nii.gz"
    )
    input_mni = "%s%s%s%s" % (
        FSLDIR,
        "/data/standard/MNI152_T1_",
        vox_size,
        "_brain.nii.gz",
    )
    input_mni_mask = "%s%s%s%s" % (
        FSLDIR,
        "/data/standard/MNI152_T1_",
        vox_size,
        "_brain_mask.nii.gz",
    )
    input_mni_sched = "%s%s" % (FSLDIR, "/etc/flirtsch/T1_2_MNI152_2mm.cnf")

    t1w_brain = output_path / f"sub-{subject}_ses-{ses}_T1w_brain.nii.gz"
    t1w_normalized = output_path / f"sub-{subject}_ses-{ses}_T1w_normalized.nii.gz"
    t1w_brain_aligned = (
        output_path / f"sub-{subject}_ses-{ses}_T1w_brain_aligned.nii.gz"
    )
    t12mni_xfm_init = output_path / f"sub-{subject}_ses-{ses}_t12mni_xfm_init.mat"
    t12mni_xfm = output_path / f"sub-{subject}_ses-{ses}_t12mni_xfm.mat"
    warp_t1w2mni = output_path / f"sub-{subject}_ses-{ses}_warp-t1w2mni.mat"

    # if not t1w_brain.exists():

    # Normalize
    print("\nRunning Normalization")
    mgru.normalize_t1w(input_t1w, t1w_normalized)

    # Skull stripping
    print("\nRunning 3dSkullStrip")
    mgru.t1w_skullstrip(t1w_normalized, str(t1w_brain))

    # Voxel reshape
    print("\nRunning Voxel Matching")
    match_target_vox_res(str(t1w_brain))

    # Create linear transform/ initializer T1w-->MNI
    print("\nInitial T1w->MNI transform")
    mgru.align(
        t1w_brain,
        input_mni,
        xfm=t12mni_xfm_init,
        bins=None,
        interp="spline",
        out=None,
        dof=12,
        cost="mutualinfo",
        searchrad=True,
    )

    # Registration from t1w -> MNI
    if nonlinear:
        print("\nRunning non-linear registration: T1w-->MNI ...")
        # Use FNIRT to nonlinearly align T1 to MNI template
        mgru.align_nonlinear(
            t1w_brain,
            input_mni,
            xfm=t12mni_xfm_init,
            out=t1w_brain_aligned,
            warp=warp_t1w2mni,
            ref_mask=input_mni_mask,
            config=input_mni_sched,
        )
    else:
        # Falling back to linear registration
        print("\nRunning linear registration: T1w-->MNI ...")
        mgru.align(
            t1w_brain,
            input_mni,
            xfm=t12mni_xfm,
            init=t12mni_xfm_init,
            bins=None,
            dof=12,
            cost="mutualinfo",
            searchrad=True,
            interp="spline",
            out=t1w_brain_aligned,
            sch=None,
        )

    # Segment wm, gm, csf
    print("\nSegmenting brain regions")
    maps = mgru.segment_t1w(t1w_brain_aligned, output_path / "masks" / f"sub-{subject}")

    wm_mask = maps["wm_prob"]
    gm_mask = maps["gm_prob"]
    csf_mask = maps["csf_prob"]
    match_target_vox_res(wm_mask)
    match_target_vox_res(gm_mask)
    match_target_vox_res(csf_mask)


from argparse import ArgumentParser


def main():
    """Starting point of the pipeline, assuming that you are using a BIDS organized dataset
    """
    parser = ArgumentParser(
        description="This is a registration script for structural MRIs."
    )
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
        "participant_label",
        help="The label(s) of the "
        "participant(s) that should be analyzed. The label "
        "corresponds to sub-<participant_label> from the BIDS "
        'spec (so it does not include "sub-"). If this '
        "parameter is not provided all subjects should be "
        "analyzed. Multiple participants can be specified "
        "with a space separated list.",
    )
    parser.add_argument(
        "--session_label",
        help="The label(s) of the "
        "session that should be analyzed. The label "
        "corresponds to ses-<participant_label> from the BIDS "
        'spec (so it does not include "ses-"). If this '
        "parameter is not provided all sessions should be "
        "analyzed. Multiple sessions can be specified "
        "with a space separated list.",
        default=1,
    )
    parser.add_argument(
        "--vox",
        default="1mm",
        help="Voxel size : 2mm, 1mm. Voxel size to use for template registrations.gi",
    )
    parser.add_argument(
        "--nonlinear", default=False, help="Whether to use nonlinear registration"
    )

    result = parser.parse_args()

    inDir = result.input_dir
    outDir = result.output_dir
    subj = result.participant_label
    sesh = result.session_label
    nonlinear = result.nonlinear
    vox = result.vox

    register_t1w_2_mni(inDir, outDir, subj, sesh, nonlinear, vox)


if __name__ == "__main__":
    main()
