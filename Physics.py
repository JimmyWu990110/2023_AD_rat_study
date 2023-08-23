import math

import numpy as np

class Physics:
    def __init__(self):
        self.gyro_hz = 42.576 # MHz/T
        self.gyro_rad = self.gyro_hz * 2 * math.pi # 10e6 rad/(s*T)
        
    def get_resonance_freq(self, B0):
        # in MHz
        return self.gyro_hz * B0
    
    def ppm2hz(self, B0, offset):
        resonance_freq = self.get_resonance_freq(B0) # in MHz
        return resonance_freq * offset # freq in hz

    def hz2ppm(self, B0, freq):
        resonance_freq = self.get_resonance_freq(B0)
        return freq / resonance_freq # offset in ppm

