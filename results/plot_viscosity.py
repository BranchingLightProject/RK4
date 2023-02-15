import numpy as np
import matplotlib.pyplot as pl
import seaborn as sns
from scipy import stats
import scipy.optimize as opt

sns.set_theme(style='whitegrid')

# old nu
data = """
1.092   4.3
1.201   8.175
1.330  9.362
1.480  10.62
1.661  11.475
"""

data = """
1.116   4.3
1.263   8.175
1.441  9.362
1.661  10.62
1.984  11.475
"""

gl = []
s_len = []

for i, line in enumerate(data.strip().split('\n')):
    line = line.split()
    gl.append(float(line[0]))
    s_len.append(float(line[1]))

y_errors_linear = [
#3.15175444627274,
#0.783575126710898,
#2.07370663750801,
#3.90230242948134,
#2.98358839869041,
1.51252606662629,
0.783575126710898,
2.07370663750801,
1.8961336357968,
1.81599076693688,

]

x_errors = [0.434*0.01/x for x in gl]
y_errors = [0.434*y_errors_linear[i]/y for i, y in enumerate(s_len)]

gl = np.log(gl)
s_len = np.log(s_len)

coef = np.polyfit(gl[:], s_len[:], 1)
f_fit = np.poly1d(coef)

coef2, V = np.polyfit(gl[1:], s_len[1:], 1, cov=True)
f_fit2 = np.poly1d(coef2)

r_value_1 = stats.linregress(gl[:], s_len[:])[2]
r_value_2 = stats.linregress(gl[1:], s_len[1:])[2]
print(r_value_1**2)
print(r_value_2**2)
print(f_fit)

op_p = opt.curve_fit(lambda x, b: 0.75*x + b, gl[1:], s_len[1:])[0]

theo = lambda x: 0.75*x + op_p[0]

fit_label = f'$m={coef2[0]:.3f}({int(np.sqrt(V[0][0])*2000)})$\n$b={coef2[1]:.3f}({int(np.sqrt(V[1][1])*2000)})$\n$r²={r_value_2**2:.3f}$'

pl.plot(gl[1:], f_fit2(gl[1:]), c='darkred', label=fit_label, linestyle="--")
pl.errorbar(gl, s_len, color='black', xerr=x_errors, yerr=y_errors, fmt='o')
pl.xlabel(r'$\log(\nu)$')

pl.ylabel(r'$\log (d_0)$')
pl.legend()

pl.subplots_adjust(left=0.1, right=0.9, top=0.9, bottom=0.1)
pl.savefig('test.png')


