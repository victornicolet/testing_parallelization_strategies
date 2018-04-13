import sys

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

data = pd.read_csv(sys.argv[1])

data_by_thread = {}

for i in [1,2,3,4,5,6,7,8,12,16,20,24,32]:
    df = data.loc[data['threads'] == i].drop(['threads'], axis=1)
    df['split/fused'] = df.apply(lambda x: x['split'] / x['fused'], axis=1)
    df = df.drop(['sequential', 'split', 'fused'], axis=1)
    df = df.pivot(index='height', columns='width',values='split/fused')
    data_by_thread[i] = df


def plot_df(df):
    heatmap = plt.pcolor(df)
    plt.title("Running time comparison for split and fused loops.")
    plt.yticks(np.arange(0.5, len(df.index), 1), df.index)
    plt.xticks(np.arange(0.5, len(df.columns), 1), df.columns)
    leg = plt.colorbar(heatmap)
    leg.ax.set_ylabel("Elapsed time split / fused, height vertical.")
    plt.show()

plot_df(data_by_thread[int(sys.argv[2])])