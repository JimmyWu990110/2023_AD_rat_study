
from Initialization import Initialization
from Lorentz_fitting import Lorentz_fitting
from Visualization import Visualization
from EMR_fitting import EMR_fitting

def fit_Lorentz_2pool_onepixel(offset, Zspec):
    initialization = Initialization(B0_shift=0)
    initialization.init_2_pool()
    fitting = Lorentz_fitting(initialization)
    # fit the data 
    fitted_paras = fitting.fit_2_pool(offset, Zspec)
    y_estimated = fitting.generate_Zpsec("2_pool", offset, fitted_paras)
    print("offset:", offset)
    print("Zspec:", Zspec)
    # print("fitted parameters:", fitted_paras)
    visualization = Visualization()
    visualization.format_paras(fitted_paras)
    visualization.plot_2_Zspec(offset, Zspec, y_estimated, ["real", "fitted"],
                               3.5)
    visualization.plot_component_2pool(offset, fitted_paras)
    visualization.plot_component_2pool(offset[8:], fitted_paras)

def fit_Lorentz_5pool_onepixel(offset, Zspec):
    initialization = Initialization(B0_shift=0)
    initialization.init_5_pool()
    fitting = Lorentz_fitting(initialization)
    # fit the data 
    fitted_paras = fitting.fit_5_pool(offset, Zspec)
    y_estimated = fitting.generate_Zpsec("5_pool", offset, fitted_paras)
    print("offset:", offset)
    print("Zspec:", Zspec)
    # print("fitted parameters:", fitted_paras)
    visualization = Visualization()
    visualization.format_paras(fitted_paras)
    visualization.plot_2_Zspec(offset, Zspec, y_estimated, ["real", "fitted"],
                               3.5)
    visualization.plot_component_5pool(offset, fitted_paras)
    visualization.plot_component_5pool(offset[8:], fitted_paras)
    
def fit_EMR_onepixel(offset, Zspec):
    freq = offset * 500
    model_type = 0
    T1w_obs = 2
    T2w_obs = 0.04
    fitting = EMR_fitting(T1w_obs=T1w_obs, T2w_obs=T2w_obs)
    fitting.set_model_type(model_type)
    initialization = Initialization(B0_shift=0)
    initialization.init_EMR(T1w_obs/T2w_obs, model_type)
    fitting.set_x0(initialization.x0)
    fitting.set_lb(initialization.lb)
    fitting.set_ub(initialization.ub)
    # fit the data 
    fitted_paras, y_estimated = fitting.fit(freq[:8], Zspec[:8], constrained=True)
    # fitted_paras, y_estimated = fitting.fit(freq, Zspec, constrained=True)
    all_paras = fitting.cal_paras()
    print("All parameters:\n", all_paras)
    visualization = Visualization()
    # visualization.plot_2_Zspec(offset, Zspec, y_estimated, ["real", "fitted"],
    #                            3.5)
    Zspec_all_extrap = fitting.generate_Zpsec(freq, fitted_paras)
    visualization.plot_2_Zspec(offset, Zspec, Zspec_all_extrap, ["real", "fitted"]
                               )





    
    
    
    