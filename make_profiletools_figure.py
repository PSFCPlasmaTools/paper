from __future__ import division

import matplotlib
matplotlib.rc('font', size=20)

import scipy
import eqtools
import profiletools
import gptools
import matplotlib.pyplot as plt
plt.ion()
plt.close('all')

shot = 1101014006
t_min = 0.965
t_max = 1.365

e = eqtools.CModEFITTree(shot)
p_CTS = profiletools.TeCTS(shot, t_min=t_min, t_max=t_max, efit_tree=e)
p_CTS.time_average(weighted=True)
p_ETS = profiletools.TeETS(shot, t_min=t_min, t_max=t_max, efit_tree=e)
p_ETS.time_average(weighted=True)
p_GPC = profiletools.TeGPC(shot, t_min=t_min, t_max=t_max, efit_tree=e)
p_GPC.time_average()
p_GPC2 = profiletools.TeGPC2(shot, t_min=t_min, t_max=t_max, efit_tree=e)
p_GPC2.time_average()

f, sl = e.plotFlux(fill=False)
sl.set_val(47)
f.axes[0].plot(p_CTS.X[:, 0], p_CTS.X[:, 1], 'gs', markersize=12, label='CTS')
f.axes[0].plot(p_ETS.X[:, 0], p_ETS.X[:, 1], 'rs', markersize=12, label='ETS')
f.axes[0].plot(p_GPC2.X[:, 0], scipy.zeros_like(p_GPC2.y), 'm^', markersize=12, label='GPC2')
f.axes[0].plot(p_GPC.X[:, 0], scipy.zeros_like(p_GPC.y), 'b^', markersize=12, label='GPC')
f.canvas.draw()

p_CTS.convert_abscissa('r/a')
p_ETS.convert_abscissa('r/a')
p_GPC.convert_abscissa('r/a')
p_GPC.remove_edge_points()
p_GPC2.convert_abscissa('r/a')
p_GPC2.remove_edge_points()

fp = plt.figure()
ap = fp.add_subplot(1, 1, 1)
p_CTS.plot_data(ax=ap, color='g', marker='s', ls='', label='CTS')
p_ETS.plot_data(ax=ap, color='r', marker='s', ls='', label='ETS')
p_GPC.plot_data(ax=ap, color='b', marker='^', ls='', label='GPC')
p_GPC2.plot_data(ax=ap, color='m', marker='^', ls='', label='GPC2')

p = p_CTS

p.add_profile(p_ETS)
p.add_profile(p_GPC)
p.add_profile(p_GPC2)

p.create_gp()

b = 10
m = 1
be = 10
me = 0.5
p.gp.k.hyperprior = (
    gptools.UniformJointPrior([(0, 30)]) *
    gptools.GammaJointPrior(
        [1 + m * b, 1 + me * be, 1, 1 + 1.01 * 200],
        [b, be, 1/0.1, 200]
    )
)

p.find_gp_MAP_estimate(verbose=True)
roa = scipy.linspace(0, 1.2)
p.smooth(roa, plot=True, ax=ap, label='smoothed profile')

ap.set_xlim(left=0, right=1.1)
ap.set_ylim(bottom=0, top=3.5)
ap.set_ylabel('$T_e$ [keV]')
ap.legend(loc='best')
fp.canvas.draw()
