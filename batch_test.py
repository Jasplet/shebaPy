#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 24 16:52:19 2021

@author: ja17375

Trail batch processing. 
"""
from batch import serial_process
path= '/Users/ja17375/Projects/Matisse_Synthetics/ppv1/ideal'
with open(f'{path}/Inversion_synthetics_noise_0.01.events', 'r+') as infile:
    events = [line.strip() + '.BH*' for line in infile]
    
test = events[0:2]
result = serial_process(test, f'{path}/run', ['SKS','SKS'])
print(result)