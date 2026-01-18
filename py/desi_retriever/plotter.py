import matplotlib.pyplot as plt
import scipy.stats
import numpy as np


def plot(D, model=None, fig=None, percs=None, title=None, merged=False):
    if fig is None:
        fig = plt.gcf()
        plt.clf()
    axes = []
    if not merged:
        for i in range(3):
            axes.append(fig.add_subplot(311 + i))
    else:
        axes = plt.gca()
    colHash = {'b': 'blue', 'r': 'green', 'z': 'red'}
    colHash1 = {'b': 'lightblue', 'r': 'lightgreen', 'z': 'orange'}
    for i, filt in enumerate('brz'):
        spec = D[filt + '_flux'] * 1
        spec[spec == 0] = np.nan
        if not merged:
            cur_ax = axes[i]
        else:
            cur_ax = axes
        cur_ax.plot(D[filt + '_wavelength'], spec, color=colHash[filt])
        if title is not None and i == 0:
            cur_ax.set_title(title)

    if model is not None:
        for i, filt in enumerate('brz'):
            spec = model[filt + '_model'] * 1
            spec[spec == 0] = np.nan
            if not merged:
                cur_ax = axes[i]
            else:
                cur_ax = axes
            cur_ax.plot(model[filt + '_wavelength'],
                        spec,
                        color=colHash1[filt],
                        alpha=0.5)
    if percs is not None:
        min_vals = []
        max_vals = []
        for i, filt in enumerate('brz'):
            vals = [
                scipy.stats.scoreatpercentile(D[filt + '_flux'], _)
                for _ in percs
            ]
            min_vals.append(vals[0])
            max_vals.append(vals[1])

            if not merged:
                axes[i].set_ylim(*vals)
        if merged:
            axes.set_ylim(min(min_vals), max(max_vals))
    plt.xlabel('Wavelength')
    plt.draw()
