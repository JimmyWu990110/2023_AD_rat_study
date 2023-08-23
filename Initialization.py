import numpy as np

class Initialization:
    def __init__(self, B0_shift): 
        self.B0_shift = B0_shift
        self.x0 = None
        self.ub = None 
        self.ub = None
    
    def init_2_pool(self):
        """
        [amp_mt, pos_mt, width_mt, 
        amp_ds, pos_ds, width_ds]
        amp is not in percentile, but in [0, 1]!
        """
        x0 = [40/100, -2.34, 25,
              40/100, 0, 2]
        lb = [0, -3, 0,
              0, -0.2, 0]
        ub = [100/100, -1.5, 9999,
              100/100, 0.2, 10]
        for i in range(2):
            x0[3*i + 1] += self.B0_shift
            lb[3*i + 1] += self.B0_shift
            ub[3*i + 1] += self.B0_shift
        self.x0 = np.array(x0)
        self.lb = np.array(lb)
        self.ub = np.array(ub)
            
    def init_5_pool(self):
        """
        [amp_noe, pos_noe, width_noe,
         amp_mt, pos_mt, width_mt, 
         amp_ds, pos_ds, width_ds,
         amp_cest, poscest, width_cest, 
         amp_apt, pos_apt, width_apt]
        amp is not in percentile, but in [0, 1]!
        """
        x0 = [5/100, -3.5, 1,
              40/100, -2.34, 25,
              40/100, 0, 2,
              3/100, 2, 1,
              5/100, 3.5, 1]
        lb = [0, -3.75, 0,
              0, -3, 0,
              0, -0.2, 0,
              0, 1.75, 0,
              0, 3.25, 0]
        ub = [10/100, -3.25, 10,
              100/100, -1.5, 9999,
              100/100, 0.2, 10,
              10/100, 2.25, 10,
              10/100, 3.75, 10]
        for i in range(5):
            x0[3*i + 1] += self.B0_shift
            lb[3*i + 1] += self.B0_shift
            ub[3*i + 1] += self.B0_shift
        self.x0 = np.array(x0)
        self.lb = np.array(lb)
        self.ub = np.array(ub)
    
    def init_EMR_model_0(self, T1w_obs_T2w_obs_ratio):
        # [R, R*M0m*T1w, T1w/T2w, T2m]
        self.x0 = [20, 2, T1w_obs_T2w_obs_ratio, 10e-6]
        self.lb = [0, 0, 0, 0]
        self.ub = [50, 10, 50, 100e-6]
        
    def init_EMR_model_1(self, T1w_obs_T2w_obs_ratio):
        # [R, M0m*T1w, T1w/T2w, T2m]
        self.x0 = [20, 0.1, T1w_obs_T2w_obs_ratio, 10e-6]
        self.lb = [0, 0, 0, 0]
        self.ub = [50, 0.5, 50, 100e-6]
            
    def init_EMR_model_2(self, T1w_obs_T2w_obs_ratio):
        # [M0m*T1w, T1w/T2w, T2m]
        self.x0 = [0.1, T1w_obs_T2w_obs_ratio, 10e-6]
        self.lb = [0, 0, 0]
        self.ub = [0.5, 50, 100e-6]
        
    def init_EMR(self, T1w_obs_T2w_obs_ratio, model_type):
        if model_type == 0:
            self.init_EMR_model_0(T1w_obs_T2w_obs_ratio)
        elif model_type == 1:
            self.init_EMR_model_1(T1w_obs_T2w_obs_ratio)
        elif model_type == 2:
            self.init_EMR_model_2(T1w_obs_T2w_obs_ratio)
    
    
    














