import os

from plot import generate_plot, heatmap, scatter

ROOT_DIR = '.'

BEAM_DIR = f"{ROOT_DIR}/beams/"
POT_DIR = f"{ROOT_DIR}/potentials/"
CORR_DIR_1D = f"{ROOT_DIR}/correlation/"
CORR_DIR_2D = f"{ROOT_DIR}/corr2d/"
SCINT_DIR = f"{ROOT_DIR}/scint/"

OPTIONS = {
    '1':{
        'dir':BEAM_DIR,
        'f':heatmap,
        'kwargs':{'x_label':'x [mm]', 'y_label':'y [mm]', 'cb_label':'|P(x,y)|', 'laser':True}
    },
    '2':{
        'dir':POT_DIR,
        'f':heatmap,
        'kwargs':{'x_label':'x [mm]', 'y_label':'y [mm]', 'cb_label':'V(x,y)'}
    },
    '3':{
        'dir':CORR_DIR_1D,
        'f':scatter,
        'kwargs':{'x_label':'dr [mm]', 'y_label':'c(dr)'}
    },
    '4':{
        'dir':SCINT_DIR,
        'f':scatter,
        'kwargs':{'x_label':'x [mm]', 'y_label':'S(x)'}
    },
    '5':{
        'dir':CORR_DIR_2D,
        'f':heatmap,
        'kwargs':{'x_label':'dx [mm]', 'y_label':'dy [mm]', 'cb_label':'c(dx, dy)', 'corr':True}
    }
}

msg = "Enter the dict you want plotted.\nOptions are:\n1: beams\n2: potentials\n3. correlation function 1D\n4. scintillation index\n5. correlation function 2D\n"

selection = OPTIONS[input(msg)]

dir_to_plot = selection['dir']
files = [f for f in os.listdir(dir_to_plot) if f.split('.')[-1] == 'csv']


for file in files:
    generate_plot(dir_to_plot+file, selection['f'], **selection['kwargs'])

    print(f"Plotting {file}...")

    os.system("gnuplot plot.gnu")

    os.system("rm plot.gnu")
