from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import sys

n = sys.argv[1]
propiedad = sys.argv[2]


plt.rcParams['figure.figsize'] = (7,5)

with open('../propiedades') as f:
    T = f.read()

def cubic(x,a,b,c,d):
    return a*x**3+b*x**2+c*x+d

df = pd.read_csv('gergdata')
compare_df = pd.read_csv('compare').dropna()
conc_df = pd.read_csv('tmp_conc')
conc_df = conc_df[conc_df['conc'] != 0].set_index('comp')
conc_df['conc'] = conc_df['conc'].round(4)*100
concentraciones = conc_df.to_string()

print('Concentraciones')
print(concentraciones)

X = compare_df.iloc[:,0].values
Y = compare_df.iloc[:,1].values

params = curve_fit(cubic,X,Y)[0]

x = df[(df['P'] > X.min()) & (df['P'] < X.max())]['P'].values
y = cubic(x,*params)

io = df[df.P == x[0]].index[0]
yhat = df[propiedad][io:io+len(x)]

error = np.abs(yhat-y)/y
error = error.mean()*100

print('Datos de comparación')
print(compare_df)
print('Datos GERG')
print(df.loc[yhat.index])

print(f'Incertidumbre: {error}')

fig, ax = plt.subplots()

ax.plot(x,y, label='Experimental (ajuste)', color='gray', ls='--',alpha=0.8)
ax.scatter(X,Y, label='Experimental', marker='x', color='black')
ax.plot(x,yhat, label=f'GERG (incertidumbre: {round(error,2)}%)')
ax.set_title(f'Velocidad del sonido {concentraciones} a T={T}K')
ax.set_xlabel('P (MPa)')
ax.set_ylabel('w (m/s)')
ax.legend()

ax.annotate(concentraciones, (200,150), xycoords='figure pixels')
plt.tight_layout()

plt.savefig(f'figs/{n}.png',dpi=200)