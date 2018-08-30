import pybimodal
import scipy.io


def convert_values_to_branch(mat_dict):
    mat_dict.copy()
    for k, v in mat_dict.items():
        try:
            expanded = [v[i, 0] for i in range(v.shape[0])]
        except IndexError:
            expanded = v
        mat_dict[k] = [pybimodal.Branch(br) for br in expanded]
    return mat_dict


mat_fbcont = scipy.io.loadmat('amplitude_continuation_emcs_pruned.mat')
fbamp_branches = {'phase0': mat_fbcont['famp_cur560_phase0_pruned'], 'phasePi2': mat_fbcont['famp_cur560_phasePi2_pruned']}
fbamp_branches = convert_values_to_branch(fbamp_branches)

mat_phcont = scipy.io.loadmat('phase_continuations.mat')
phases = mat_phcont['rephases'].T
phase_branches = {
    '0.0': phases[:3],
    '0.2': phases[3:6],
    '0.34': phases[6:8],
    '0.36': phases[8],
    '0.38': phases[9]
}
phase_branches = convert_values_to_branch(phase_branches)
