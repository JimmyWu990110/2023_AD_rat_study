import math

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from Data import DataPreprocessing, DataLoader
from Physics import Physics
from Visualization import Visualization
from APTwCalculation import APTwCalculation


def test_DataPreprocessing():
    base_dir = r"C:\Users\jwu191\Desktop\Projects\rat fitting"
    # case_name = "20221003_172029_AD_1_16"
    case_name = "20221003_184311_AD_1_17"
    preprocessing = DataPreprocessing(base_dir, case_name)
    preprocessing.get_scan_info()
    preprocessing.to_nifti()
    # preprocessing.read_T1(seq_id=5) # T2: 4, T1:5
    # preprocessing.read_T2(seq_id=4) # T2: 4, T1:5
    
def test_DataPreprocessing_2(): # T1, T2 map
    base_dir = r"C:\Users\jwu191\Desktop\Projects\rat fitting"
    # case_name = "20221003_172029_AD_1_16"
    case_name = "20221003_184311_AD_1_17"
    preprocessing = DataPreprocessing(base_dir, case_name)
    preprocessing.get_scan_info()
    preprocessing.read_T1_map(5)

def test_DataLoader(): # single point
    x = 41
    y = 22
    base_dir = r"C:\Users\jwu191\Desktop\Projects\rat fitting"
    case_name = "20221003_172029_AD_1_16"
    # case_name = "20221003_184311_AD_1_17"
    dataloader = DataLoader(base_dir, case_name)
    img_1 = dataloader.read_APT(group="A", B1=0.7)
    img_2 = dataloader.read_APT(group="A", B1=1.3)
    img_3 = dataloader.read_APT(group="A", B1=2)
    img_4 = dataloader.read_APT(group="A", B1=4)
    y1 = dataloader.get_Zspec_onepixel(img_1, x, y)
    offset, Zspec_1 = dataloader.preprocess_Zspec_onepixel(y1)
    y2 = dataloader.get_Zspec_onepixel(img_2, x, y)
    offset, Zspec_2 = dataloader.preprocess_Zspec_onepixel(y2)
    y3 = dataloader.get_Zspec_onepixel(img_3, x, y)
    offset, Zspec_3 = dataloader.preprocess_Zspec_onepixel(y3)
    y4 = dataloader.get_Zspec_onepixel(img_4, x, y)
    offset, Zspec_4 = dataloader.preprocess_Zspec_onepixel(y4)
    visualization = Visualization()
    visualization.plot_multi_Zspec(offset, [Zspec_1, Zspec_2, Zspec_3, Zspec_4],
                                   labels=["0.7uT", "1.3uT", "2uT", "4uT"])
    visualization.plot_multi_Zspec(offset[8:], [Zspec_1[8:], Zspec_2[8:], Zspec_3[8:], Zspec_4[8:]],
                                   labels=["0.7uT", "1.3uT", "2uT", "4uT"])
    

def test_DataLoader_2(): # roi by mask
    base_dir = r"C:\Users\jwu191\Desktop\Projects\rat fitting"
    # case_name = "20221003_172029_AD_1_16"
    case_name = "20221003_184311_AD_1_17"
    dataloader = DataLoader(base_dir, case_name)
    img_1 = dataloader.read_APT(group="A", B1=0.7)
    img_2 = dataloader.read_APT(group="A", B1=1.3)
    img_3 = dataloader.read_APT(group="A", B1=2)
    img_4 = dataloader.read_APT(group="A", B1=4)
    roi = dataloader.get_ROI_by_mask()
    y1 = dataloader.get_avg_Zspec_by_roi(img_1, roi)
    offset, Zspec_1 = dataloader.preprocess_Zspec_onepixel(y1)
    y2 = dataloader.get_avg_Zspec_by_roi(img_2, roi)
    offset, Zspec_2 = dataloader.preprocess_Zspec_onepixel(y2)
    y3 = dataloader.get_avg_Zspec_by_roi(img_3, roi)
    offset, Zspec_3 = dataloader.preprocess_Zspec_onepixel(y3)
    y4 = dataloader.get_avg_Zspec_by_roi(img_4, roi)
    offset, Zspec_4 = dataloader.preprocess_Zspec_onepixel(y4)
    print(100 * Zspec_2)
    visualization = Visualization()
    mask = dataloader.read_mask()
    visualization.show_img_2D(mask)
    visualization.plot_multi_Zspec(offset, [Zspec_1, Zspec_2, Zspec_3, Zspec_4],
                                   labels=["0.7uT", "1.3uT", "2uT", "4uT"])
    visualization.plot_multi_Zspec(offset[8:], [Zspec_1[8:], Zspec_2[8:], Zspec_3[8:], Zspec_4[8:]],
                                   labels=["0.7uT", "1.3uT", "2uT", "4uT"])


def test_Physics():
    physics = Physics()
    print(physics.get_resonance_freq(B0=3))
    print(physics.get_resonance_freq(B0=7))
    print(physics.get_resonance_freq(B0=11.7))
    
    a = physics.ppm2hz(B0=3, offset=3.5)
    b = physics.hz2ppm(B0=3, freq=a)
    print(a, b)

def test_APTwCalculation():
    x = 41
    y = 22
    base_dir = r"C:\Users\jwu191\Desktop\Projects\rat fitting"
    # case_name = "20221003_172029_AD_1_16"
    case_name = "20221003_184311_AD_1_17"
    dataloader = DataLoader(base_dir, case_name)
    emr = dataloader.read_APT(group="A", B1=0.7)
    wassr = dataloader.read_WASSR(group="A")
    calculation = APTwCalculation(emr, wassr, threshold=10)
    # calculation.cal_B0_shift_onepixel(wassr[:, x, y])
    B0_shift_map = calculation.cal_B0_shift_map()
    dataloader.save_arr_as_nifti(B0_shift_map, r"C:\Users\jwu191\Desktop\B0_shift_map_A.nii")
    APTw = calculation.cal_APTw()
    dataloader.save_arr_as_nifti(APTw, r"C:\Users\jwu191\Desktop\APTw_A_0p7uT.nii")



# test_DataPreprocessing()
# test_DataPreprocessing_2()
# test_DataLoader()
test_DataLoader_2()
# test_Physics()
# test_APTwCalculation()


