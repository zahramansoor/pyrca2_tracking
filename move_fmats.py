# -*- coding: utf-8 -*-
"""
Created on Mon Jan 30 10:59:23 2023

@author: Han
"""

import os, shutil
src = r"Z:\sstcre_imaging\e201"
dst = r"Z:\sstcre_analysis"
# move all converted fmats to separate folder
for i in os.listdir(src):
    pth = os.path.join(src, str(i))
    if "week" not in i: #week folders don't have sbx file
        imgfl = [os.path.join(pth, xx) for xx in os.listdir(pth) if "000" in xx][0]
    else:
        imgfl = pth
    mat = os.path.join(imgfl, "suite2p", "plane0", "Fall.mat") 
    if os.path.exists(mat):
        try:
            shutil.copy(mat, os.path.join(dst, f"e200_day{int(i):03d}_Fall.mat"))
        except: # for weeks
            shutil.copy(mat, os.path.join(dst, f"e200_week{int(i[4:]):02d}_Fall.mat"))