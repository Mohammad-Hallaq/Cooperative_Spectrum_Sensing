# -*- coding: utf-8 -*-
"""
Created on Sun Jul  4 02:11:26 2021

@author: mohammad.hq
"""
import numpy as np
from gnuradio import gr
from numpy import linalg as LA
import scipy as sp
from scipy import linalg as sLA
class blk (gr.sync_block):  # other base classes are basic_block, decim_block, interp_block
    def __init__ (self, N = 1000, L=8, inp = [0]):  # only default arguments here
          gr.sync_block.__init__ (
            self,
            name='MME Block',  
            in_sig=[np.complex64],
            out_sig=[np.float32]
             )
        self.N = N
        self.L = L
        self.inp = inp
    def work (self, input_items, output_items):
        while len(self.inp) < self.N and all(input_items[0]) != 0:
            self.inp.extend(input_items[0])
            output_items[0][:] = 0
            return len(output_items[0])
        if len(self.inp) >= self.N:
            Y = self.inp[0:self.N]
            Y_pad = np.pad(Y, (0,self.L), 'constant')
            Y_mid = Y_pad.reshape(self.L+self.N,1)
            Y_mid = Y_mid[0,:]
            auto = np.correlate(Y_mid, Y_mid, mode='full')
            auto = auto[auto.size//2:]/self.N
            R = sLA.toeplitz(auto[0:self.L])
            w,v = LA.eig(R)
            Test = max(w)/min(w)            
            output_items[0][:] = Test.real
            self.inp = []
            return len(output_items[0])
        else:
            output_items[0][:] = 0
            return len(output_items[0])
