import nipype.interfaces.fsl as fsl


def nn_downsample(in_file, out_file, xfm_file, isoxfm=1.25, verbose=False):
    nifti = ".nii.gz"
    mat = ".mat"

    if not in_file.endswith(nifti):
        in_file += nifti
    if not out_file.endswith(nifti):
        out_file += nifti
    if not xfm_file.endswith(mat):
        xfm_file += mat

    flt = fsl.FLIRT(
        in_file=in_file,
        reference=in_file,
        # out_file=out_file,
        out_matrix_file=xfm_file,
        interp="nearestneighbour",
        apply_isoxfm=isoxfm,
    )
    flt_res = flt.run()

    # dilm = fsl.maths.DilateImage(
    #     in_file=flt_res.outputs.out_file,
    #     operation="mean",
    #     kernel_shape="3D",
    #     out_file=out_file,
    # )
    # dilm.run()
