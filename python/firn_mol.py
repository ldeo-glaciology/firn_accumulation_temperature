# Main solver functions_classes

def setup(self,**kwargs):
        # define defaults
        self.params['T'] = 10
        self.params['r_s_sim'] = 'default'

        self.params['b0_mpy'] = 0.1
        self.params['beta'] = 1
        self.params['nu'] = 1
        self.params['T_s_dim'] = 253.15
        self.params['plotting'] = 1
        self.params['saving_xt'] = 1
        self.params['plotting_period'] = 3000
        self.params['sampling_period'] = 20
        self.params['dz'] = 0.01
        self.params['t_total'] = 20
        self.params['save_dir'] = 'results_scratch'
        self.params['save'] = 1
        self.params['sim_T'] = True
        self.params['sim_r'] = True
        self.params['PauseGrainEvolution_t'] = np.nan
        self.params['z0'] = 100
        self.params['r_s_dim'] = 2.5000e-07   # grain size at the surface (0.5 mm)^2
        self.params['phi_s'] = 0.5
        self.params['n'] = 1
        self.params['simDuration'] = 10
        self.params['scaleDuration'] = 0




    #import numpy as np
    #p  = {'n': 1.0, 'beta': np.arange(0,10,1)}



    print(simDuration)


    return b0_mpy

def can_you_run_scripts_from_same_file():
    import firn_mol
    firn_mol.setup()


def run(p):
    n = p['n']
    beta = p['beta']
    model_results = n*beta
    return model_results
