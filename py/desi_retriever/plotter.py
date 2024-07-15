import matplotlib.pyplot as plt
import scipy.stats
import numpy as np


def plot(D, model=None, fig=None, percs=None, title=None):
    if fig is None:
        fig = plt.gcf()
        plt.clf()
    axes = []
    for i in range(3):
        axes.append(fig.add_subplot(311 + i))
    colHash = {'b': 'blue', 'r': 'green', 'z': 'red'}
    colHash1 = {'b': 'lightblue', 'r': 'lightgreen', 'z': 'orange'}
    for i, filt in enumerate('brz'):
        spec = D[filt + '_flux'] * 1
        spec[spec == 0] = np.nan
        axes[i].plot(D[filt + '_wavelength'],
                     D[filt + '_flux'],
                     color=colHash[filt])
        if title is not None and i == 0:
            axes[i].set_title(title)

    if model is not None:
        for i, filt in enumerate('brz'):
            axes[i].plot(model[filt + '_wavelength'],
                         model[filt + '_model'],
                         color=colHash1[filt],
                         alpha=0.5)
    if percs is not None:
        for i, filt in enumerate('brz'):
            vals = [
                scipy.stats.scoreatpercentile(D[filt + '_flux'], _)
                for _ in percs
            ]
            axes[i].set_ylim(*vals)
    plt.xlabel('Wavelength')
    plt.draw()
