import os 
from os import path
import shutil
import glob

dirs=next(os.walk('.'))[1]
cwd=path.abspath("./")
case_dirs=[]
for folder in dirs:
    subdir=path.join(cwd,folder)
    if "processor0" in next(os.walk(subdir))[1]:
        case_dirs.append(subdir)

for case_dir in case_dirs:
    processors = glob.iglob(os.path.join(case_dir, "processor*"))
    for processor in processors:
        if path.isdir(processor):
            shutil.rmtree(processor)
