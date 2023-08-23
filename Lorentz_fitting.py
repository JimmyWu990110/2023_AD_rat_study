import os

import numpy as np
from scipy.optimize import curve_fit


class Lorentz_fitting:
    def __init__(self, initialization): 
        self.initialization = initialization
    
    def lorentz(self, x, amp, pos, width):
        factor = ((x - pos) / width) * ((x - pos) / width)
        return amp / (1 + 4*factor)
    
    def func_2_pool(self, x, amp_mt, pos_mt, width_mt,
                    amp_ds, pos_ds, width_ds):
        mt = self.lorentz(x, amp_mt, pos_mt, width_mt)
        ds = self.lorentz(x, amp_ds, pos_ds, width_ds)
        return 1 - (mt + ds)
    
    def func_5_pool(self, x, amp_noe, pos_noe, width_noe,
                   amp_mt, pos_mt, width_mt,
                   amp_ds, pos_ds, width_ds, 
                   amp_cest, pos_cest, width_cest,
                   amp_apt, pos_apt, width_apt):
        noe = self.lorentz(x, amp_noe, pos_noe, width_noe)
        mt = self.lorentz(x, amp_mt, pos_mt, width_mt)
        ds = self.lorentz(x, amp_ds, pos_ds, width_ds)
        cest = self.lorentz(x, amp_cest, pos_cest, width_cest)
        apt = self.lorentz(x, amp_apt, pos_apt, width_apt)
        return 1 - (noe + mt + ds + cest + apt)
    
    def fit_2_pool(self, offset, Zspec):
        """
        2 pool model (MT, DS), given offset and Zspec,
        return fitted parameters (shape: [6,])
        """
        popt, pcov = curve_fit(self.func_2_pool, xdata=offset, ydata=Zspec, 
                               p0=self.initialization.x0,
                               bounds=(self.initialization.lb, self.initialization.ub), 
                               method="trf", maxfev=5000)
        for i in range(2):
            popt[3*i + 1] -= self.initialization.B0_shift
        return popt

    def fit_5_pool(self, offset, Zspec):
        """
        5 pool model (nOE, MT, DS, CEST@2ppm, APT), given offset and Zspec,
        return fitted parameters (shape: [15,])
        """    
        popt, pcov = curve_fit(self.func_5_pool, xdata=offset, ydata=Zspec, 
                               p0=self.initialization.x0,
                               bounds=(self.initialization.lb, self.initialization.ub), 
                               method="trf", maxfev=5000)
        for i in range(5):
            popt[3*i + 1] -= self.initialization.B0_shift
        return popt
    
    def generate_Zpsec(self, model_type, offset, paras):
        """
        given offset and fitted parameters, calculate and return the Zspec
        model type should be given
        return Zspec (shape: [n,]), where n is the num of offsets
        """
        y = []
        if model_type == "2_pool":
            for x in offset: 
                y.append(self.func_2_pool(x, *paras))
        elif model_type == "5_pool":
            for x in offset: 
                y.append(self.func_5_pool(x, *paras))
        else:
            raise Exception("Invalid model type!")
        return np.array(y)
    
        
    
    
    def func_fixed_5_pool(self, x, amp_noe, width_noe, amp_mt, width_mt, 
                         amp_ds, width_ds, amp_cest, width_cest, amp_apt, width_apt):
        noe = self.lorentz(x, amp_noe, -3.5, width_noe)
        mt = self.lorentz(x, amp_mt, -2.34, width_mt)
        ds = self.lorentz(x, amp_ds, 0, width_ds)
        cest = self.lorentz(x, amp_cest, 2, width_cest)
        apt = self.lorentz(x, amp_apt, 3.5, width_apt)
        return 100 - (noe + mt + ds + cest + apt)
    



    

  










