import os
from dmriprep.utils.bids import get_bids_layout
from dmriprep.workflows.dwi.base import init_dwi_preproc_wf
from dmriprep.workflows.dwi.util import init_dwi_concat_wf

from joblib import Parallel, delayed
from pathlib import Path

def run_dmriprep(sub):
    bdir = "/mnt/ssd3/j1c/hcp1200"
    out_dir = "/mnt/ssd3/j1c/hcp1200_processed"

    ses = 1
    sub_dict = get_bids_layout(bdir, sub, ses)
    if len(sub_dict[ses].keys()) == 1:
        dwi_file = sub_dict[ses][1]["dwi_file"]
        fbvec = sub_dict[ses][1]["fbvec"]
        fbval = sub_dict[ses][1]["fbval"]
        metadata = sub_dict[ses][1]["metadata"]
    else:
        dwi_files = []
        fbvecs = []
        fbvals = []
        metadata_files = []
        for acq in sub_dict[ses].keys():
            dwi_files.append(sub_dict[ses][acq]["dwi_file"])
            fbvecs.append(sub_dict[ses][acq]["fbvec"])
            fbvals.append(sub_dict[ses][acq]["fbval"])
            metadata_files.append(sub_dict[ses][acq]["metadata"])
        wf = init_dwi_concat_wf(dwi_files, fbvals, fbvecs, sub, ses, out_dir)
        out = wf.run()
    wf = init_dwi_preproc_wf(sub, ses, dwi_file, fbval, fbvec, metadata, out_dir, vox_size="1mm")
    out = wf.run()
    
    
    
path = Path("../hcp1200")
p = [x.parts[-1][4:] for x in path.iterdir()]

out = Parallel(n_jobs=20, verbose=1)(delayed(run_dmriprep)(s) for s in p[:20])