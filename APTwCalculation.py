import os

import numpy as np
import matplotlib.pyplot as plt
from scipy import interpolate

class APTwCalculation:
    
    def __init__(self, threshold=0, min_=-15, max_=0):
        self.threshold = threshold
        self.min = min_
        self.max = max_
        # 43 offsets
        self.offset = np.array([200, 200, 
                                4, -4, 3.75, -3.75, 3.5, -3.5, 3.5, -3.5, 3.5, -3.5, 
                                3.5, -3.5, 3.5, -3.5, 3.5, -3.5, 3.25, -3.25, 3, -3, 
                                2.5, -2.5, 2, -2, 1.5, -1.5, 1, -1, 0.5, -0.5, 0.25, -0.25, 
                                0, 6, 8, 10, 12, 14, 16, 18, 20])
        # 22 offsets
        self.WASSR_offset = np.array([40, 
                                      1, 0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 
                                      0.1, 0, -0.1, -0.2, -0.3, -0.4, -0.5, -0.6, 
                                      -0.7, -0.8, -0.9, -1])
    
    def cal_zero_mask(self, EMR):
        M0 = EMR[1]
        zero_mask = np.ones((M0.shape[0], M0.shape[1]))
        for i in range(M0.shape[0]):
            for j in range(M0.shape[1]):
                if M0[i][j] < self.threshold:
                    zero_mask[i][j] = 0
        return zero_mask
    
    def cal_B0_shift_onepixel(self, WASSR_Spec):
        # poly fit the WASSR Spec and find the min point of the upsampled offset
        sort_index = np.argsort(self.WASSR_offset)
        x_sorted = self.WASSR_offset[sort_index] * 500 # to hz (for 11.7 T)
        y_sorted = WASSR_Spec[sort_index]
        m0 = y_sorted[-1]
        # normalized by M0, but fitting not include M0
        paras = np.polyfit(x_sorted[:-1], y_sorted[:-1]/m0, deg=12)
        p = np.poly1d(paras)
        x_upsampled = np.arange(-500, 501, 1)
        index = np.argmin(p(x_upsampled))
        # plt.scatter(x_sorted[:-1], y_sorted[:-1]/m0)
        # plt.plot(x_upsampled, p(x_upsampled))
        # plt.show()
        return x_upsampled[index]
                
    def cal_B0_shift_map(self, WASSR, zero_mask):
        B0_shift_map = np.zeros((WASSR.shape[1], WASSR.shape[2]))
        for i in range(WASSR.shape[1]):
            for j in range(WASSR.shape[2]):
                if zero_mask[i][j] != 0:
                    B0_shift_map[i][j] = self.cal_B0_shift_onepixel(WASSR[:, i, j])
        return B0_shift_map
    
    def cal_MTRsaym_onepixel(self, Zspec, B0_shift):
        # print("B0_shift:", B0_shift)
        m0 = Zspec[1]
        Zspec /= m0 # normalization
        freq = self.offset * 500 # to hz (for 11.7 T)
        # print(freq)
        # print(Zspec)
        pos_1500_hz = np.mean(Zspec[np.where(freq == 1500)[0]]) # 3 ppm
        pos_1750_hz = np.mean(Zspec[np.where(freq == 1750)[0]]) # 3.5 ppm
        pos_2000_hz = np.mean(Zspec[np.where(freq == 2000)[0]]) # 4 ppm
        neg_1500_hz = np.mean(Zspec[np.where(freq == -1500)[0]])
        neg_1750_hz = np.mean(Zspec[np.where(freq == -1750)[0]])
        neg_2000_hz = np.mean(Zspec[np.where(freq == -2000)[0]])
        # calculate positive 3.5ppm(corrected)
        x_pos = np.array([1500, 1750, 2000])
        y_pos = np.array([pos_1500_hz, pos_1750_hz, pos_2000_hz])
        x_interp_pos = np.arange(1750-500, 1750+500+1, 1) # max shift: -/+ 500 hz
        func_pos = interpolate.interp1d(x_pos, y_pos, "linear", fill_value="extrapolate")
        y_interp_pos = func_pos(x_interp_pos)
        pos_1750_hz_corrected = y_interp_pos[np.where(x_interp_pos == 1750+B0_shift)][0]
        # plt.plot(x_interp_pos, y_interp_pos)
        # plt.scatter(x_pos, y_pos)
        # plt.show()  
        # calculate negative 3.5ppm(corrected)
        x_neg = np.array([-2000, -1750, -1500])
        y_neg = np.array([neg_2000_hz, neg_1750_hz, neg_1500_hz])
        x_interp_neg = np.arange(-1750-500, -1750+500+1, 1) # max shift: -/+ 500 hz
        func_neg = interpolate.interp1d(x_neg, y_neg, "linear", fill_value="extrapolate")
        y_interp_neg = func_neg(x_interp_neg)
        neg_1750_hz_corrected = y_interp_neg[np.where(x_interp_neg == -1750+B0_shift)][0]
        # plt.plot(x_interp_neg, y_interp_neg)
        # plt.scatter(x_neg, y_neg)
        # plt.show()   
        # print("MTR asym before correction:", 100 * (neg_1750_hz-pos_1750_hz))
        # print("MTR asym after correction:", 100 * (neg_1750_hz_corrected-pos_1750_hz_corrected))
        return 100 * (neg_1750_hz_corrected - pos_1750_hz_corrected)
        # return 100 * (neg_1750_hz-pos_1750_hz)
               
    def cal_APTw(self, WASSR, EMR):
        """
        calculate the APTw img using the given EMR and WASSR
        EMR shape: [43 64, 64], WASSR shape: [22, 64, 64]
        first calculate the zero mask using EMR
        then calculate the B0_shift_map using WASSR
        finally calculate APTw using EMR and B0_shift_map
        return B0_shift_map (shape: [64, 64]) and APTw (shape: [64, 64])
        """
        print("********", "cal_APTw", "********")
        zero_mask = self.cal_zero_mask(EMR)
        B0_shift_map = self.cal_B0_shift_map(WASSR, zero_mask)
        APTw = np.zeros((EMR.shape[1], EMR.shape[2]))
        for i in range(EMR.shape[1]):
            for j in range(EMR.shape[2]):
                if zero_mask[i][j] == 0 or EMR[:, i, j][1] == 0:
                    APTw[i][j] = self.min
                else:
                    APTw[i][j] = self.cal_MTRsaym_onepixel(EMR[:, i, j], B0_shift_map[i][j])
        APTw[np.where(APTw < self.min)] = self.min
        APTw[np.where(APTw > self.max)] = self.max
        return B0_shift_map, APTw
    
        



    


