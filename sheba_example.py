#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 18 16:45:38 2021

@author: ja17375

Simple example (for one phase) of how to use the SHEBA wrapper
"""
import obspy
from wrapper import Wrapper

data = obspy.read('example/data/HKT_1995226_043721_SKS.BH?')

Example = Wrapper(data, 'SKS', rundir='example/')
# Preprocess method filters and trims data
Example.preprocess()
result = Example.measure_splitting('HKT_example', sheba_exec_path='/Users/ja17375/Ext_programs/bin', window=False)
result.to_csv('HKT_example.sdb', index=False, sep=' ')