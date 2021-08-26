#!/usr/bin/env python3

import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import sys
import os


def plot_ij(prop):
    for ij in df['ij'].unique():
        plt.gcf().clf()
        msk = (df['ij'] == ij)
        mean_u = df[msk].u.mean()
        print(ij, prop, "mean_u: ", mean_u)

        if mean_u > 0.01:
            print('[CHECK]',prop,ij,mean_u)
            # Comparar puntos según como x,y
            x = df[msk].rho
            yhat = df[msk][prop]
            y = df[msk][f'test_{prop}']
            # Para hacer diagrama de linealidad
            #x = df[msk][f'test_{prop}']
            #y = df[msk][prop]
            plt.scatter(x, yhat, label=ij)
            plt.scatter(x, y, color='gray', alpha=0.8, marker='x')

            plt.title(rf'${prop}$ según $\rho$')
            plt.xlabel(r'Densidad $\rho$ (mol/L)')
            plt.ylabel(f'${prop}$')
            plt.grid()
            plt.legend()
            plt.savefig(f'figs/{prop}/{ij}.png')
        else:
            print('[OK]',prop,ij,mean_u)

file = sys.argv[1]

df = pd.read_csv(file, sep=' ')

variables = ['P','cv', 'w', 'cp']

for var in variables:
    df['u'] = 100*np.abs(df[f'test_{var}']-df[var])/df[var]
    df['ij'] = df['i'].astype(str) + ',' + df['j'].astype(str)
    df[f'delta_{var}'] = np.abs(df[var] - df[f'test_{var}'])
    os.system(f'echo Removing old figs...;rm -v ./figs/{var}/*')

    plot_ij(var)

df.to_csv('deltas.csv', sep='\t')
