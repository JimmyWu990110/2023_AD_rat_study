import os

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import SimpleITK as sitk
from brukerapi.dataset import Dataset
from scipy.io import loadmat


class DataPreprocessing:
    
    def __init__(self, base_dir, case_name):
        self.base_dir = base_dir
        self.case_name = case_name
        self.output_dir = self._create_output_dir()
        
    def _create_output_dir(self):
        os.makedirs(os.path.join(self.base_dir, "processed_data", self.case_name), 
                    exist_ok=True)
        return os.path.join(self.base_dir, "processed_data", self.case_name)
                
    def _sort_seq_id(self, seq_ids):
        sorted_ids = [] # list of str -> list of int
        for seq_id in seq_ids:
            sorted_ids.append(int(seq_id))
        sort_index = np.argsort(sorted_ids)
        return sort_index
    
    def _add_seq_detail(self, data):
        seq_ids = data[:, 0]
        protocals = data[:, 1]
        APT = ["APT_A_1p3uT", "APT_A_0p7uT", "APT_A_2uT", "APT_A_4uT",
               "APT_B_1p3uT"]
        WASSR = ["WASSR_A", "WASSR_B"]
        APT_counter = 0
        WASSR_counter = 0
        for i in range(data.shape[0]):
            if "_APT" in protocals[i]:
                data[i, 1] = APT[APT_counter]
                APT_counter += 1
            if "_WASSR" in protocals[i]:
                data[i, 1] = WASSR[WASSR_counter]
                WASSR_counter += 1
        return data
        
    def get_scan_info(self):
        """
        read folder_detail (given), sort the items by seq_id, add seq details,
        and finally write the results in "all_scan_info.xlsx"
        """
        print("********", "get_scan_info", "********")
        df = pd.read_excel(os.path.join(self.base_dir, "raw_data", self.case_name, 
                                        "folder_detail.xlsx"))
        data = np.row_stack(([np.array(df.columns)], np.array(df)[:-1])) # [N, 2]
        seq_ids = data[:, 0]
        sort_index = self._sort_seq_id(seq_ids)
        data = data[sort_index]
        data = self._add_seq_detail(data)
        print("scan info:\n", data)
        column_names = ["Sequence_id", "Protocol"]
        df_result = pd.DataFrame(data=data, columns=column_names)
        df_result.to_excel(os.path.join(self.output_dir, "all_scan_info.xlsx"), 
                           index=False)
    
    def _reshape_img_3D(self, img):
        # eg. for APT: [64, 64, N] -> [N, 64, 64]
        img_reshaped = np.zeros((img.shape[2], img.shape[0], img.shape[0]))
        for i in range(img_reshaped.shape[0]):
            img_reshaped[i] = img[:, :, i]
        return img_reshaped
    
    def _rotate_img_3D(self, img):
        tmp = np.rot90(img, k=1) # k: num of rot90 (counterclockwise)
        return np.flip(tmp, axis=1) # axis=0: up-down, axis=1: left-right    
    
    def read_bruker_img(self, seq_id, pid=1):
        """
        sequence id, pdata id (1 or 2)
        """
        dataset = Dataset(os.path.join(self.base_dir, "raw_data", self.case_name, 
                                       str(seq_id), "pdata", str(pid)))
        img = np.squeeze(dataset.data)
        # first rotate, then reshape! (for [64, 64, N], it's axes (0->1) by default)
        img = self._rotate_img_3D(img)
        return self._reshape_img_3D(img)
    
    def _normalize_T1_map(self, T1_map):
        T1_map[T1_map < 1185] = 1185
        T1_map[T1_map > 3815] = 3815
        return T1_map
    
    def _normalize_T2_map(self, T2_map):
        T2_map[T2_map < 0] = 0
        T2_map[T2_map > 80] = 80
        return T2_map
        
    def read_T1_map(self, seq_id):
        dataset = Dataset(os.path.join(self.base_dir, "raw_data", self.case_name, 
                                       str(seq_id), "pdata", "2"))
        img = np.squeeze(dataset.data) # [64, 64, 5, 11]
        for i in range(5):
            tmp = img[:, :, i, :]
            tmp = self._rotate_img_3D(tmp)
            self._save_arr_as_nifti(self._reshape_img_3D(tmp), 
                                    os.path.join(self.output_dir, "1_nifti", 
                                                 "T1_map_"+str(i)+".nii"))
        data = loadmat(os.path.join(self.base_dir, "raw_data", self.case_name, 
                                    "T1map.mat"))
        # img = data["T1_Map"]
        # img = self._rotate_img_3D(img)
        # return self._reshape_img_3D(img)  
        T1_map = self._rotate_img_3D(img[:, :, 2, :])
        T1_map = self._reshape_img_3D(T1_map) 
        return self._normalize_T1_map(T1_map)
 

    def read_T2_map(self, seq_id):
        dataset = Dataset(os.path.join(self.base_dir, "raw_data", self.case_name, 
                                       str(seq_id), "pdata", "2"))
        img = np.squeeze(dataset.data) # [64, 64, 5, 11]
        for i in range(5):
            tmp = img[:, :, i, :]
            tmp = self._rotate_img_3D(tmp)
            self._save_arr_as_nifti(self._reshape_img_3D(tmp), 
                                    os.path.join(self.output_dir, "1_nifti", 
                                                 "T2_map_"+str(i)+".nii"))  
        data = loadmat(os.path.join(self.base_dir, "raw_data", self.case_name, 
                                    "T2map.mat"))
        # img = data["T2_Map"]
        # img = self._rotate_img_3D(img)
        # return self._reshape_img_3D(img) 
        T2_map = self._rotate_img_3D(img[:, :, 2, :])
        T2_map = self._reshape_img_3D(T2_map)  
        return self._normalize_T2_map(T2_map)         

    def _save_arr_as_nifti(self, arr, path):
        """
        given an arr, save it to given path in nifti format
        """
        # print("save to", path, "shape:", arr.shape)
        img = sitk.GetImageFromArray(arr)
        sitk.WriteImage(img, path)
    
    def to_nifti(self):
        """
        read sequences info from "all_scan_info.xlsx", 
        save APT, WASSR and T1, T2 maps in nifti files under "1_nifti" folder
        """
        print("********", "to_nifti", "********")
        os.makedirs(os.path.join(self.output_dir, "1_nifti"), exist_ok=True)
        df = pd.read_excel(os.path.join(self.output_dir, "all_scan_info.xlsx"))
        data = np.array(df)
        for i in range(data.shape[0]):
            if "APT_" in data[i, 1] or "WASSR_" in data[i, 1]:
                img = self.read_bruker_img(seq_id=data[i, 0]) # type: int
                filename = str(data[i, 0]) + "_" + data[i, 1] + ".nii" # id_protocal
                path = os.path.join(self.output_dir, "1_nifti", filename)
                self._save_arr_as_nifti(img, path)
            if data[i, 1].startswith("T1map"):
                T1_map = self.read_T1_map(seq_id=data[i, 0])
                self._save_arr_as_nifti(T1_map, os.path.join(self.output_dir, "1_nifti", 
                                                             "T1_map.nii"))            
            if data[i, 1].startswith("T2map"):
                T2_map = self.read_T2_map(seq_id=data[i, 0])
                self._save_arr_as_nifti(T2_map, os.path.join(self.output_dir, "1_nifti", 
                                                             "T2_map.nii"))
                

class DataLoader:
    
    def __init__(self, base_dir, case_name):
        self.base_dir = base_dir
        self.case_name = case_name
        self.nifti_dir = os.path.join(self.base_dir, "processed_data", 
                                      self.case_name, "1_nifti")
        # APT: 43 offsets, 32 unique
        self.offset = np.array([200, 200, 
                                4, -4, 3.75, -3.75, 3.5, -3.5, 3.5, -3.5, 3.5, -3.5, 
                                3.5, -3.5, 3.5, -3.5, 3.5, -3.5, 3.25, -3.25, 3, -3, 
                                2.5, -2.5, 2, -2, 1.5, -1.5, 1, -1, 0.5, -0.5, 0.25, -0.25, 
                                0, 6, 8, 10, 12, 14, 16, 18, 20])
        self.unique_offset = np.array([200, 
                                       4, -4, 3.75, -3.75, 3.5, -3.5, 3.25, -3.25, 
                                       3, -3, 2.5, -2.5, 2, -2, 1.5, -1.5, 1, -1, 
                                       0.5, -0.5, 0.25, -0.25, 
                                       0, 6, 8, 10, 12, 14, 16, 18, 20])
        # WASSR: 22 offsets
        self.WASSR_offset = np.array([40, 
                                      1, 0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 
                                      0.1, 0, -0.1, -0.2, -0.3, -0.4, -0.5, -0.6, 
                                      -0.7, -0.8, -0.9, -1])
        
    def read_APT(self, group, B1):
        """
        read APT image given group (A or B) and B1 power, 
        return arr (not img) (shape: [43, 64, 64])
        """
        assert(group == "A" or group == "B"), "group can only be A or B"
        assert(B1 == 0.7 or B1 == 1.3 or B1 == 2 or B1 == 4), "B1 = 0.6, 1.3, 2 or 4 uT"
        suffix = "0p7uT"
        if B1 == 1.3:
            suffix = "1p3uT"
        elif B1 == 2:
            suffix = "2uT"
        elif B1 == 4:
            suffix = "4uT"
        for f in os.listdir(self.nifti_dir):
            if group+"_"+suffix in f:               
                img = sitk.ReadImage(os.path.join(self.nifti_dir, f))
                # print(os.path.join(self.nifti_dir, f))
                return sitk.GetArrayFromImage(img)
        raise Exception("Data not found")
        
    def read_WASSR(self, group):
        """
        read WASSR image given group (A or B), return arr (not img)
        return the WASSR array (shape: [22, 64, 64])
        """
        assert(group == "A" or group == "B"), "group can only be A or B"
        for f in os.listdir(self.nifti_dir):
            if "WASSR_"+group in f:               
                img = sitk.ReadImage(os.path.join(self.nifti_dir, f))
                # print(os.path.join(self.nifti_dir, f))
                return sitk.GetArrayFromImage(img)
        raise Exception("Data not found")
        
    def read_mask(self):
        """
        read the ROI mask, (2D, 64*64, only one of the 43 dynamics)
        get the slice (dynamic first)
        return the mask array (shape: [64, 64])
        """
        filename = self.case_name + "_mask.nii"
        mask = sitk.ReadImage(os.path.join(self.nifti_dir, filename), sitk.sitkInt32)
        return sitk.GetArrayFromImage(mask)[0]
        
    def get_ROI_by_mask(self, label=9999):
        """
        If given a label (nonzero), return the ROI corresponding to that label,
        If not given a label, return all ROIs (all area with nonzero label)
        return roi (shape:)
        """
        filename = self.case_name + "_mask.nii"
        mask = sitk.ReadImage(os.path.join(self.nifti_dir, filename), sitk.sitkInt32)
        mask_arr = self.read_mask()
        print("mask shape:", mask_arr.shape)
        roi = np.where(mask_arr != 0)
        if label != 9999:
            roi = np.where(mask_arr == label)
        print("mask size:", roi[0].shape[0])
        return roi
        
    def get_Zspec_onepixel(self, img, x, y):
        """
        given coordinates x and y, get the raw Zspec from the [43, 64, 64] img,
        return raw Zspec (shape: [43,])
        """
        return img[:, x, y]
    
    def get_avg_Zspec_by_roi(self, img, roi):
        """
        given an image and corresponding roi, 
        return the average raw Zspec (shape: [43,])
        """
        Zspec = []
        for i in range(roi[0].shape[0]):
            Zspec.append(img[:, roi[0][i], roi[1][i]])
        return np.mean(np.array(Zspec), axis=0)
    
    def preprocess_Zspec_onepixel(self, y):
        """
        given a raw Zspec in onepixel (shape: [43,]), calculate avg for duplicate offsets,
        normalized by M0 (only the second one) and sort the Zspec
        return sorted Zspec without M0 as offset (shape: [32,])
        """
        Zspec = np.zeros(self.unique_offset.shape[0])
        Zspec[0] = y[1] # only use the second M0!
        Zspec[5] = np.mean(y[np.where(self.offset == 3.5)])
        Zspec[6] = np.mean(y[np.where(self.offset == -3.5)])
        for i in range(4):
            Zspec[i+1] = y[i+2]
        for i in range(25):
            Zspec[i+7] = y[i+18] 
        Zspec /= Zspec[0] # normalize by M0
        # sort offset and Zspec, from downfield to upfield
        sort_index = np.argsort(self.unique_offset)[::-1]
        return self.unique_offset[sort_index][1:], Zspec[sort_index][1:] # not include M0
        
    def _save_arr_as_nifti(self, arr, path):
        """
        given an arr, save it to given path in nifti format
        """
        # print("save to", path, "shape:", arr.shape)
        img = sitk.GetImageFromArray(arr)
        sitk.WriteImage(img, path)
    

        
    








    
    
    
        
        