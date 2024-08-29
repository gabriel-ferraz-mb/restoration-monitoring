# -*- coding: utf-8 -*-
"""
Created on Sun Jun  2 19:35:30 2024

@author: Gabriel
"""

import os
os.chdir(r'C:\Projetos\bioflore\ecosia\PSdata')

l=os.listdir()
l = [val for val in l if not val.endswith("QA.tif")]

for f in l:
    qa = f.replace('.tif', '_QA.tif')
    temp = f.replace('.tif', '_temp.tif')
    
    s_data = os.path.getsize(f)
    s_qa = os.path.getsize(qa)
    if s_qa > s_data:
        os.rename(f, temp)
        os.rename(qa, f)
        os.rename(temp, qa)
        