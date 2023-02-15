import numpy as np
import matplotlib.pyplot as pl
import seaborn as sns
from PIL import Image


sns.set_theme(style="darkgrid")

FOLDER = "./scint"
FILE = "heatmap_p17"
X_REAL = 36.1
SCINT_MAX_X = 0.0


def intensities(folder, file):
    raw_data = np.loadtxt(f"{folder}/{file}_reversed.csv", delimiter=",")
    _ = pl.hist(raw_data[:, 2], bins="auto")
    pl.savefig(f"{folder}/{file}_hist.png")
    pl.close()


def scintillation(folder, file, x_real=10.0, dis=0.0):
    scint = np.loadtxt(f"{folder}/{file}_plane.csv", delimiter=",")
    x_max = len(scint)
    scint_max = (0.0, 0.0)

    for i in range(x_max):
        s = scint[i]
        if s[1] > scint_max[1] and s[0] < 1000.6:
            scint_max = s

    scint_clean = np.zeros((x_max - 10,2))
    # scint_max = (0.0, 0.0)

    # for j in range(x_max - 10):
    #     i = j + 2
    #     clean_val = 0.5 * scint[i][1] + 0.1 * (
    #         scint[i - 5][1]
    #         + scint[i - 4][1]
    #         + scint[i - 3][1]
    #         + scint[i - 2][1]
    #         + scint[i - 1][1]
    #         + scint[i + 1][1]
    #         + scint[i + 2][1]
    #         + scint[i + 3][1]
    #         + scint[i + 4][1]
    #         + scint[i + 5][1]
    #     )
    #     if clean_val > scint_max[1] and scint[i][0]<0.6:
    #         scint_max = (scint[i][0], clean_val)
    #     scint_clean[j][0] = scint[i][0]
    #     scint_clean[j][1] = clean_val

    x = np.linspace(0, x_real, num=(x_max - 10))
    scint_y = np.linspace(0.0, 1.1 * scint_max[1])
    scint_x = np.ones_like(scint_y)

    print("scint max in: ", scint_max[0])
    print(scint_max[1])

    pl.plot(scint[:,0], scint[:,1])
    # pl.plot(scint_clean[:,0], scint_clean[:,1])
    pl.plot(scint_max[0] * scint_x, scint_y, "r", linewidth=2)

    # pl.text(x[scint_max[0]], scint_max[1], f"({x[scint_max[0]]}, {scint_max[1]})")
    pl.text(
        1.05 * scint_max[0], 1.01 * scint_max[1], f"$d_0 = {scint_max[0]:.2f}$ mm"
    )
    pl.xlabel("x (mm)")
    pl.ylabel("S")
    pl.tight_layout(pad=0.5)
    pl.savefig(f"{folder}/{file}_scint.png")


def merge_images(folder, file_clean):
    scint = Image.open(f"{folder}/{file_clean}_scint.svg")
    branches = Image.open(f"{folder}/{file_clean}.png")

    scint = scint.resize(branches.size)

    merged = Image.new("RGB", (branches.size[0], 2 * branches.size[1]), (250, 250, 250))
    merged.paste(branches, (0, 0))
    merged.paste(scint, (0, branches.size[1]))
    merged.save(f"{folder}/{file_clean}_merged.png", "PNG")


def main():
    # intensities(FOLDER, FILE)
    scintillation(FOLDER, FILE, x_real=10.0)
    # merge_images(FOLDER, FILE)


if __name__ == "__main__":
    main()
