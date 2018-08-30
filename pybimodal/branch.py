import numpy as np
import pandas as pd


def variable_names():
    return ['Re(Es)', 'Im(Es)', '|Ew|^2', 'rho', 'n']


def param_names():
    return ['kappa_s', 'kappa_w', 'mu_s', 'mu_w', 'epsi_ss', 'epsi_ww', 'epsi_sw', 'epsi_ws', 'beta', 'J_p', 'eta',
            'tau_r', 'S_in', 'V', 'Z_QD', 'n_bg', 'tau_sp', 'T_2', 'A', 'hbar_omega', 'epsi_tilda', 'J', 'feed_phase',
            'feed_ampli', 'tau_fb', 'epsi0', 'hbar', 'e0', 'alpha_par', 'omega1']


def param_plot_names():
    return ['\\kappa_s', '\\kappa_w', '\\mu_s', '\\mu_w', '\\epsilon_{ss}', '\\epsilon_{ww}', '\\epsilon_{sw}',
            '\\epsilon_{ws}', '\\beta', 'J_p', '\\eta', '\\tau_r', 'S^{in}', 'V', 'Z^{QD}', 'n_{bg}', '\\tau_{sp}',
            'T_2', 'A', 'hbar\\omega', '\\epsilon_0n_{bg}c_{0}', 'J', 'Feedback Phase', 'Feedback Amp', '\\tau_{fb}',
            '\\epislon0', 'hbar', 'e_0', '\\alpha', '\\omega1', '\\omega2']


class Branch:
    def __init__(self, matlab_branch):
        points = matlab_branch['point'][0, 0]
        x = np.concatenate(points['x'].squeeze(), axis=1).astype(np.float).T
        params = np.concatenate(points['parameter'].squeeze(), axis=0).astype(np.float)
        nunst = matlab_branch['nunst'][0, 0]
        ind_hopf = matlab_branch['indHopf'][0, 0] - 1
        ind_fold = matlab_branch['indFold'][0, 0] - 1

        hopf, fold = np.zeros_like(nunst, dtype=np.bool), np.zeros_like(nunst, dtype=np.bool)
        hopf[ind_hopf, :] = True
        fold[ind_fold, :] = True
        intensity_x1x2 = np.square(np.linalg.norm(x[:, 1:2], axis=1))[:, np.newaxis]

        self.variable_names = variable_names()
        self.param_names = param_names()
        self.extended_names = ['|Es|^2', 'nunst', 'hopf', 'fold']

        data = np.concatenate([x, params, intensity_x1x2, nunst, hopf, fold], axis=1)
        col_names = self.variable_names + self.param_names + self.extended_names
        self._df = pd.DataFrame(data=data, columns=col_names)

        self._active_params = self._get_active_params()

    def _get_active_params(self):
        names = []
        for name in self.param_names:
            if np.allclose(self._df.loc[0, name], self._df.loc[:, name]):
                pass
            else:
                names.append(name)
        return self.df.columns.intersection(names)

    @property
    def df(self):
        return self._df

    @property
    def cont_params(self):
        return self._active_params

    @property
    def cont(self):
        return self._df.loc[:, self.cont_params]

    @property
    def vars(self):
        return self._df.loc[:, self.variable_names]

    @property
    def x(self):
        cols = self.cont_params.drop('omega1', errors='ignore')
        cols = cols.insert(loc=0, item='|Es|^2')
        return self._df.loc[:, cols]
