import matplotlib.pyplot as plt
import pandas as pd
import sys


def compare(i, j, df):
    plt.gcf().clf()
    fig, ax = plt.subplots()
    ax2 = ax.twinx()
    df_ij = ij(i, j, df)
    df_ij[['tau2', 'ar', 'ao']].plot(ax=ax)
    df_ij[['cv']].plot(ax=ax2, color='red')
    ax2.scatter(x=df_ij.index, y=df_ij['cv_test'], label='test_data', color='red')

    plt.legend()

    return df_ij


def ij(i, j, df):
    return df[(df['i'] == i) & (df['j'] == j)]


df = pd.read_csv(sys.argv[1])
R = 8.3144720

df['cv'] = -df['tau']**2*(df['ao']+df['ar'])
df[['cv']] = df[['cv']]*R
df['ao'] = -df['ao']
df['ar'] = -df['ar']

for i in range(1, 22):
    for j in range(1, 22):
        compare(i, j, df)
        plt.savefig(f'figs/{i},{j}.png')
