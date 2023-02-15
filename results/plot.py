"""
Generate gnuplot scripts
"""

def generate_plot(filename, function, **kwargs):
    filename = '.'.join( filename.split('.')[:-1] )

    string = f"""set datafile separator ','\nset term png\nset output '{filename}.png'\n{function(filename, **kwargs)}"""

    with open("plot.gnu", 'w+') as f:
        f.write(string)

def heatmap(filename, x_label='', y_label='', cb_label='', ratio=False, laser=False, corr=False):
    if ratio: r = "set size ratio 1"
    else: r = ''
    if laser: laser_str = """set palette defined (0 'black', 0.3 'green' , 1 'white')
unset colorbox"""
    else: laser_str = ''
    
    if corr: range_ = """set xrange[-1:1]
set yrange [-1:1]"""
    else: range_ = """set xrange[0:10]
set yrange [0:5]"""


    string = f"""set pm3d map
{r}
{range_}
set xlabel '{x_label}'
set ylabel '{y_label}'
set cblabel '{cb_label}'
{laser_str}

splot '{filename}.csv' notitle"""

    return string

def vector2D(filename):
    string = f"plot '{filename}.csv' w vec"

    return string

def scatter(filename, x_label='', y_label=''):
    string = f"""set xlabel '{x_label}'
set ylabel '{y_label}'

plot '{filename}.csv' notitle w l"""

    return string

def scatter3D(filename):
    string = f"splot '{filename}.csv' notitle"

    return string

def heat_vector(filename, ratio=False):
    if ratio: r = "set size ratio 1"
    else: r = ''

    string = f"""set pm3d map
{r}

splot '{filename}.csv' using 1:2:5 notitle, '' using 1:2:(0):3:4:(0) notitle w vec linecolor rgb '#000000'"""

    return string
