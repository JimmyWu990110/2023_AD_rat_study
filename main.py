import os

import numpy as np

from Data import DataPreprocessing, DataLoader
from Visualization import Visualization
from APTwCalculation import APTwCalculation

    
base_dir = r"C:\Users\jwu191\Desktop\Projects\rat fitting"
# case_list = ["20221003_184311_AD_1_17"]
case_list = os.listdir(os.path.join(base_dir, "raw_data"))


# B0 = 11.7 # T
# B1 = 1.3 # uT

def preprocess(base_dir, case_name):
    preprocessing = DataPreprocessing(base_dir, case_name)
    preprocessing.get_scan_info()
    preprocessing.to_nifti()

def cal_APTw(base_dir, case_name):
    output_dir = os.path.join(base_dir, "processed_data", case_name, "2_calculated")
    os.makedirs(output_dir, exist_ok=True)
    dataloader = DataLoader(base_dir, case_name)
    wassr_a = dataloader.read_WASSR(group="A")
    emr_a_0p7 = dataloader.read_APT(group="A", B1=0.7)
    emr_a_1p3 = dataloader.read_APT(group="A", B1=1.3)
    emr_a_2 = dataloader.read_APT(group="A", B1=2)
    emr_a_4 = dataloader.read_APT(group="A", B1=4)
    wassr_b = dataloader.read_WASSR(group="B")
    emr_b_1p3 = dataloader.read_APT(group="B", B1=1.3)
    calculation = APTwCalculation(threshold=10)
    B0_shift_map_A, APTw_A_0p7 = calculation.cal_APTw(wassr_a, emr_a_0p7)
    B0_shift_map_A, APTw_A_1p3 = calculation.cal_APTw(wassr_a, emr_a_1p3)
    B0_shift_map_A, APTw_A_2 = calculation.cal_APTw(wassr_a, emr_a_2)
    B0_shift_map_A, APTw_A_4 = calculation.cal_APTw(wassr_a, emr_a_4)
    B0_shift_map_B, APTw_B_1p3 = calculation.cal_APTw(wassr_b, emr_b_1p3)
    dataloader._save_arr_as_nifti(B0_shift_map_A, os.path.join(output_dir, "B0_shift_map_A.nii"))
    dataloader._save_arr_as_nifti(B0_shift_map_B, os.path.join(output_dir, "B0_shift_map_B.nii"))
    dataloader._save_arr_as_nifti(APTw_A_0p7, os.path.join(output_dir, "APTw_A_0p7uT.nii"))
    dataloader._save_arr_as_nifti(APTw_A_1p3, os.path.join(output_dir, "APTw_A_1p3uT.nii"))
    dataloader._save_arr_as_nifti(APTw_A_2, os.path.join(output_dir, "APTw_A_2uT.nii"))
    dataloader._save_arr_as_nifti(APTw_A_4, os.path.join(output_dir, "APTw_A_4uT.nii"))
    dataloader._save_arr_as_nifti(APTw_B_1p3, os.path.join(output_dir, "APTw_B_1p3uT.nii"))
    
    
visualization = Visualization()
# DataPreprocessing
for case_name in case_list:
    preprocess(base_dir, case_name)
    cal_APTw(base_dir, case_name)
    





# fit_Lorentz_2pool_onepixel(offset, Zspec)
# fit_Lorentz_5pool_onepixel(offset, Zspec)
# fit_EMR_onepixel(offset, Zspec)
