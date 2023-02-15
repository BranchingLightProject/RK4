import os

import numpy as np
import matplotlib.pyplot as pl
import seaborn as sns
from PIL import Image


sns.set_theme(style="darkgrid")

FOLDER = "./correlation"


def correlation(folder, file, x_real=10.0, dis=0.0):
    scint = np.loadtxt(f"{folder}/{file}.csv", delimiter=",")
    x_max = len(scint)
    pl.plot(scint[:,0], scint[:,1], label='$l_c=0.117$ mm')


    pl.xlabel("x (mm)")
    pl.ylabel("C(dr)")
    pl.tight_layout(pad=0.5)
    pl.legend()
    pl.savefig(f"{folder}/{file}_corr.png")
    pl.close()


def main():
#    files = os.listdir(FOLDER)
    files = ['heatmap_p13.csv']
    
    for file in files:
        print(file)
        if file.split('.')[-1] != 'csv': continue
        correlation(FOLDER, '.'.join(file.split('.')[:-1]), x_real=10.0)
    


if __name__ == "__main__":
    main()
