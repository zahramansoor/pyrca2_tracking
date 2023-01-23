# -*- coding: utf-8 -*-
"""
Created on Mon Jan 23 10:47:19 2023

@author: Han
"""

import os, numpy as np, tifffile as tif

dr = r'H:\imaging_data'
concattifs = [] #in order of days?
days = os.listdir(dr); days.sort()
for day in :
    daydr = os.path.join(dr, daydr)
    # find images from session, split up in 3000 frame tifs
    tifs = [os.path.join(daydr, tif) for tifs in os.listdir(daydr) if '.tif' in tif]