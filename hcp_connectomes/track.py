# external package imports
import nibabel as nib
import numpy as np
from dipy.core.gradients import gradient_table
from dipy.data import get_sphere
from dipy.direction import ProbabilisticDirectionGetter, peaks_from_model
from dipy.io.gradients import read_bvals_bvecs
from dipy.reconst.csdeconv import ConstrainedSphericalDeconvModel, recursive_response
from dipy.reconst.shm import CsaOdfModel
from dipy.tracking import utils
from dipy.tracking.local_tracking import LocalTracking
from dipy.tracking.stopping_criterion import BinaryStoppingCriterion
from dipy.tracking.streamline import Streamlines


def build_seed_list(wm_mask, dens):
    """uses dipy tractography utilities in order to create a seed list for tractography
    Parameters
    ----------
    wm_mask : np.array
    dens : int
        seed density
    Returns
    -------
    ndarray
        locations for the seeds
    """
    stream_affine = np.eye(4)

    seeds = utils.random_seeds_from_mask(
        wm_mask, affine=stream_affine, seeds_count=int(dens), seed_count_per_voxel=True
    )
    return seeds


def odf_mod_est(gtab):
    print("Fitting CSA ODF model...")
    mod = CsaOdfModel(gtab, sh_order=6)
    return mod


def csd_mod_est(dwi, gtab, wm_mask):
    print("Fitting CSD model...")
    print("Estimating recursive response...")
    response = recursive_response(
        gtab,
        dwi,
        mask=wm_mask,
        sh_order=6,
        peak_thr=0.01,
        init_fa=0.08,
        init_trace=0.0021,
        iter=8,
        convergence=0.001,
        parallel=False,
    )
    mod = ConstrainedSphericalDeconvModel(gtab, response, sh_order=6)
    return mod


def load_data(fdwi, fbval, fbvec, fwmparc):
    dwi = nib.load(fdwi).get_fdata()

    bvals, bvecs = read_bvals_bvecs(fbval, fbvec)
    gtab = gradient_table(bvals, bvecs)

    wmparc = nib.load(fwmparc).get_fdata()
    # These are WM values from freesurfer
    values = list(range(251, 256)) + list(range(3000, 5003))
    wm_mask = np.isin(wmparc, values)

    return dwi, gtab, wm_mask


def run_tractography(fdwi, fbval, fbvec, fwmparc, mod_func, mod_type, seed_density=20):
    """
    mod_func : 'str'
        'csd' or 'csa'
    mod_type : 'str'
        'det' or 'prob'
    seed_density : int, default=20
        Seeding density for tractography
    """
    # Getting default params
    sphere = get_sphere("repulsion724")
    stream_affine = np.eye(4)

    # Loading data
    dwi, gtab, wm_mask = load_data(fdwi, fbval, fbvec, fwmparc)

    # Make tissue classifier
    tiss_classifier = BinaryStoppingCriterion(wm_mask)

    if mod_func == "csd":
        mod = csd_mod_est(gtab, dwi, wm_mask)
    elif mod_func == "csa":
        mod = odf_mod_est(gtab)

    # Build seed list
    seeds = utils.random_seeds_from_mask(
        wm_mask,
        affine=stream_affine,
        seeds_count=int(seed_density),
        seed_count_per_voxel=True,
    )

    # Make streamlines
    if mod_type == "det":
        print("Obtaining peaks from model...")
        mod_peaks = peaks_from_model(
            mod,
            dwi,
            sphere,
            relative_peak_threshold=0.5,
            min_separation_angle=25,
            mask=wm_mask,
            npeaks=5,
            normalize_peaks=True,
        )
        streamline_generator = LocalTracking(
            mod_peaks,
            tiss_classifier,
            seeds,
            stream_affine,
            step_size=0.5,
            return_all=True,
        )
    elif mod_type == "prob":
        print("Preparing probabilistic tracking...")
        print("Fitting model to data...")
        mod_fit = mod.fit(dwi, wm_mask)
        print("Building direction-getter...")
        try:
            print(
                "Proceeding using spherical harmonic coefficient from model estimation..."
            )
            pdg = ProbabilisticDirectionGetter.from_shcoeff(
                mod_fit.shm_coeff, max_angle=60.0, sphere=sphere
            )
        except:
            print("Proceeding using FOD PMF from model estimation...")
            fod = mod_fit.odf(sphere)
            pmf = fod.clip(min=0)
            pdg = ProbabilisticDirectionGetter.from_pmf(
                pmf, max_angle=60.0, sphere=sphere
            )
        streamline_generator = LocalTracking(
            pdg, tiss_classifier, seeds, stream_affine, step_size=0.5, return_all=True
        )

    print("Reconstructing tractogram streamlines...")
    streamlines = Streamlines(streamline_generator)
    tracks = Streamlines([track for track in streamlines if len(track) > 60])
    return tracks
