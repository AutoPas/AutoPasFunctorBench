import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

file_location: str = '/Users/irene/TUM/thesis/AutoPas_Thesis/AutoPasFunctorBench/AxilrodTeller/'

nu = 1.6152500E-3
def axilrodTeller(x:float):
    positionI = np.array([x, 0, 0])
    positionJ = np.array([1, 1, 0])
    positionK = np.array([1, -1, 0])

    R1 = positionJ - positionI
    R2 = positionK - positionJ
    R3 = positionI - positionK

    cos1 = - np.dot(R1, R3)/(np.linalg.norm(R1) * np.linalg.norm(R3))
    cos2 = - np.dot(R2, R1)/(np.linalg.norm(R2) * np.linalg.norm(R1))
    cos3 = - np.dot(R3, R2)/(np.linalg.norm(R3) * np.linalg.norm(R2))

    pot = nu * (1+3*cos1*cos2*cos3)/(np.linalg.norm(R1)*np.linalg.norm(R2)*np.linalg.norm(R3))**3

    return pot

def plot_expected_force():
    df = pd.read_csv(file_location + 'Particle1AlongX.csv', delimiter=',', index_col=False)
    ax = df.plot(x='Xposition', y='Energy')

    df['dU_dr1'] = np.gradient(df['Energy'].values, df['Xposition'].values)
    df['ForceI_x_expected'] = -df['dU_dr1'] * 1
    df['ForceI_y_expected'] = -df['dU_dr1'] * 0
    df['ForceI_z_expected'] = -df['dU_dr1'] * 0

    df.plot(x='Xposition', y='ForceI_x_expected', ax=ax, style='--')
    df.plot(x='Xposition', y='ForceI_y_expected', ax=ax, style='-.')
    df.plot(x='Xposition', y='ForceI_z_expected', ax=ax, style='--')
    plt.grid(visible=True)
    plt.xlabel('x')
    plt.legend(
        [r'$Potential$', r'$F1_{x}$', r'$F1_{y}$', r'$F1_{z}$'])
    plt.title('AT, r1=(x,0,0), r2=(1,1,0), r3=(1,-1,0) - Expected')
    #plt.savefig('Figures/31_07/AxilrodTeller/Equilateral/ExpectedForce.png')
    plt.show()

def plot_computed_force():
    df = pd.read_csv(file_location + 'Particle1AlongX.csv', delimiter=',', index_col=False)

    ax = df.plot(x='Xposition', y='Energy')
    df.plot(x='Xposition', y='Force_x', ax=ax, style='--')
    df.plot(x='Xposition', y='Force_y', ax=ax, style='-.')
    df.plot(x='Xposition', y='Force_z', ax=ax, style='--')
    plt.grid(visible=True)
    plt.xlabel('x')
    plt.legend(
        [r'$Potential$', r'$F1_{x}$', r'$F1_{y}$', r'$F1_{z}$'])
    plt.title('AT, r1=(x,0,0), r2=(1,1,0), r3=(1,-1,0) - Computed')
    #plt.savefig('Figures/31_07/AxilrodTeller/Equilateral/ComputedForce.png')
    plt.show()

def plot_potential():
    df = pd.read_csv(file_location + 'Particle1AlongX.csv', delimiter=',', index_col=False)

    ax = df.plot(x='Xposition', y='Energy', color='blue', linestyle='-')

    x: np.array = np.linspace(0.001, 1.999999999999891, 2000)
    y: np.array = [axilrodTeller(x1) for x1 in x]
    ax.plot(x, y, color='orange', linestyle='--')

    plt.legend([r'$Potential$', r'$\frac{11}{8}\frac{ \nu }{R^9}$'])
    plt.xlabel('x')
    plt.title('AT, r1=(x,0,0), r2=(1,1,0), r3=(1,-1,0)')
    plt.savefig('Figures/13_08/AxilrodTeller/Particle1AlongX/Potential.png')
    plt.show()

if __name__ == '__main__':
    plot_potential()
    plot_computed_force()
    plot_expected_force()