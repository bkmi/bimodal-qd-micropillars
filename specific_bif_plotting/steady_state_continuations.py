import pybimodal
import scipy.io
import matplotlib.pyplot as plt


def convert_values_to_branch(mat_dict):
    mat_dict.copy()
    for k, v in mat_dict.items():
        try:
            expanded = [v[i, 0] for i in range(v.shape[0])]
        except IndexError:
            expanded = v
        mat_dict[k] = [pybimodal.Branch(br) for br in expanded]
    return mat_dict


def split_nunst(branch):
    nunst_max = int(branch.df['nunst'].max())
    indices_by_nunst = {}
    for i in range(nunst_max + 1):
        df = branch.x[branch.df['nunst'] == i]
        if not df.empty:
            inds = []
            grouped_by_incremental_inds = df.groupby(df.index.to_series().diff().ne(1).cumsum())
            for v in grouped_by_incremental_inds.groups.values():
                inds.append(v)
            indices_by_nunst[i] = inds
    return indices_by_nunst


def stst_plot(branches):
    fig, ax = plt.subplots()
    for b in branches:
        inds_by_nunst = split_nunst(b)
        for k, inds in inds_by_nunst.items():
            for i in inds:
                if k == 0:
                    ax.plot(b.x.iloc[i, 1], b.x.iloc[i, 0], color='black', label='stable')
                else:
                    ax.plot(b.x.iloc[i, 1], b.x.iloc[i, 0], color='red', linestyle='dashed', label=str(k))
    return fig, ax


mat_fbcont = scipy.io.loadmat('amplitude_continuation_emcs_pruned.mat')
fbamp_branches = convert_values_to_branch(
    {'phase0': mat_fbcont['famp_cur560_phase0_pruned'], 'phasePi2': mat_fbcont['famp_cur560_phasePi2_pruned']}
)

stst_plot(fbamp_branches['phase0'])
plt.show()

mat_phcont = scipy.io.loadmat('phase_continuations.mat')
phases = mat_phcont['rephases'].T
phase_branches = convert_values_to_branch(
    {0.0: phases[:3], 0.2: phases[3:6], 0.34: phases[6:8], 0.36: phases[8], 0.38: phases[9]}
)
