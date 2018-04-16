import sys

import matplotlib.pyplot as plt
from matplotlib import cm
import numpy as np
import pandas as pd

data = pd.read_csv(sys.argv[1])
suffix=""
if len(sys.argv) > 2:
    suffix = sys.argv[2]

def plot_df(df, ncores, range_one):
    heatmap = plt.pcolor(df, cmap='RdBu', vmin=1-range_one, vmax=1 + range_one)
    plt.title("Running time comparison for split and fused loops (%d threads)" % ncores)
    plt.yticks(np.arange(0.5, len(df.index), 1), df.index)
    plt.xticks(np.arange(0.5, len(df.columns), 1), df.columns)
    leg = plt.colorbar(heatmap)
    leg.ax.set_ylabel("Elapsed time split / fused, height vertical.")
    plt.savefig("figures/mtrr_heatmap_%s_threads%s" % (ncores, suffix))
    plt.clf()



for i in [1,2,3,4,5,6,7,8,12,16,20,24,32]:
    df = data.loc[data['threads'] == i].drop(['threads'], axis=1)
    df = df.groupby(['height', 'width']).mean()
    df['split/fused'] = df.apply(lambda x: x['split'] / x['fused'], axis=1)
    df = df.drop(['sequential', 'split', 'fused'], axis=1)
    tmax = df['split/fused'].max()
    tmin = df['split/fused'].min()
    df = df.reset_index()
    df = df.pivot(index='height', columns='width',values='split/fused')
    range_one = max(abs(1-tmax), abs(1-tmin))
    plot_df(df, i, range_one)


def plot_speedups():
    df = data.loc[data['height'] == 20480].drop(['height'], axis=1)
    df.groupby(['width']).mean()
    df['speedup_split'] = df.apply(lambda x: x['sequential'] / x['split'], axis=1)
    df['speedup_fused'] = df.apply(lambda x: x['sequential'] / x['fused'], axis=1)
    df = df.drop(['sequential', 'split', 'fused'], axis=1)
    df = df.reset_index()

    for width in [4096, 8192, 12288, 16384]: #,20480,24576,28672,32768]
        plt.clf()
        dfi = df.loc[df['width'] == width].drop(['width'], axis=1)
        plt.plot(dfi['threads'], dfi['speedup_split'])
        plt.plot(dfi['threads'], dfi['speedup_fused'])
        plt.show()
        plt.title("Speedup of split and fused for different row sizes, 20480 rows x %d." % width)


# plot_speedups()


