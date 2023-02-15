import numpy as np
import matplotlib.pyplot as pl
import seaborn as sns
from scipy import stats
import scipy.optimize as opt

sns.set_theme(style='whitegrid')


data = """
4	4.3075
8	8.175
12	9.362
16	10.62
20	11.475
"""

gl = []
s_len = []

for i, line in enumerate(data.strip().split('\n')):
    line = line.split()
    gl.append(float(line[0]))
    s_len.append(float(line[1]))

y_errors_linear = [
1.51252606662629,
0.783575126710898,
2.07370663750801,
1.8961336357968,
1.81599076693688,
]

#x_errors = [0.434*0.1/x for x in gl]
#y_errors = [0.434*y_errors_linear[i]/y for i, y in enumerate(s_len)]

#gl = np.log(gl)
#s_len = np.log(s_len)

coef = np.polyfit(gl, s_len, 1)
f_fit = np.poly1d(coef)

r_value = stats.linregress(gl, s_len)[2]
print(r_value**2)
print(coef)
print(f_fit)

op_p = opt.curve_fit(lambda x, b: 0.75*x + b, gl, s_len)[0]

theo = lambda x: 0.75*x + op_p[0]

p = '('
p_ = ')'




#pl.plot(gl, f_fit(gl), c='darkred', label=f'$m={coef[0]:.3f}$\n$r²={r_value**2:.3f}$')
#pl.plot(gl, theo(gl), label='$m=0.75$')
pl.errorbar(gl, s_len, color='black', xerr=0.1, yerr=y_errors_linear, fmt='o')
pl.xlabel(r'glycerin percentage')

pl.ylabel(r'$d_0 ($mm$)$')
pl.legend()

pl.savefig('test.png')


