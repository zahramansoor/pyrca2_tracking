# -*- coding: utf-8 -*-
"""
Created on Mon Jan 30 10:59:23 2023

@author: Han
"""

import os, shutil
src = r"Z:\imaging_yc"
dst = r"Z:\cellreg1month_Fmats"
# move all converted fmats to separate folder
for i in os.listdir(src):
    pth = os.path.join(src, i)
    mat = os.path.join(pth, "suite2p", "plane0", "Fall.mat")
    if os.path.exists(mat):
        shutil.move(mat, os.path.join(dst, i+"_Fall.mat"))