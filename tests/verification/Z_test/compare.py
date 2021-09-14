from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import sys

n = sys.argv[1]

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

df.P = df.P.round(5)

X = compare_df.iloc[:,0].values
Y = compare_df.iloc[:,1].values

params = curve_fit(cubic,X,Y)[0]

x = df[(df['P'] > X.min()) & (df['P'] < X.max())]['P'].values
gerg_x = df['P']
y = cubic(x,*params)

io = df[df.P == x[0]].index[0]
yhat = df.Z[io:io+len(x)]

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
ax.plot(gerg_x, df.Z, label=f'GERG (incertidumbre: {round(error,2)}%)')
ax.set_title(f'Factor de compresibilidad {concentraciones} a T={T}K')
ax.set_xlabel('P (MPa)')
ax.set_ylabel('Z')
ax.legend()



plt.tight_layout()
plt.xlim(X.min()*0.9,X.max()*1.1)
plt.ylim(Y.min()*0.9,Y.max()*1.1)

plt.savefig(f'figs/{n}.png',dpi=200)
