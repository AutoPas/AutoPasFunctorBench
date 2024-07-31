import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# Might need to be modified
file_location: str = 'AxilrodTeller/'

class GeometryHandler:
    def __init__(self):
        self.file_name_potential: str = ''
        self.file_name_force: str = ''

        self.R1_versor: np.array = []
        self.R3_versor: np.array = []
        self.R3_versor: np.array = []

        self.dR1_dR3: float = 0

        self.title: str = ''

        self.expected_force_fig: str = ''
        self.computed_force_fig: str = ''


class EquilateralGeometryHandler(GeometryHandler):
    def __init__(self):
        super().__init__()
        self.file_name: str = 'EquilateralGeometry.csv'

        self.R1_versor: np.array = [0, 1, 0]
        self.R2_versor: np.array = [np.sqrt(3) / 2, -1 / 2, 0]
        self.R3_versor: np.array = [-np.sqrt(3) / 2, -1 / 2, 0]

        self.dR1_dR3 = 1

        self.title: str = r'Axilrod Teller $R_1 = R_2 = R_3 = R$'

        self.expected_force_fig: str = 'Figures/ExpectedForceTotal_Equilateral.png'
        self.computed_force_fig: str = 'Figures/ComputedForceTotal_Equilateral.png'

def plot_expected_force(handler: GeometryHandler):
    df = pd.read_csv(file_location + handler.file_name, delimiter=',', index_col=False)
    ax = df.plot(x='Distance', y='Energy')

    df['dU_dR1'] = np.gradient(df['Energy'].values, df['Distance'].values)
    R1_versor = handler.R1_versor
    R3_versor = handler.R3_versor
    dR1_dR3 = handler.dR1_dR3
    df['ForceI_x_expected'] = -df['dU_dR1'] * (dR1_dR3 * R3_versor[0] - R1_versor[0])
    df['ForceI_y_expected'] = -df['dU_dR1'] * (dR1_dR3 * R3_versor[1] - R1_versor[1])
    df['ForceI_z_expected'] = -df['dU_dR1'] * (dR1_dR3 * R3_versor[2] - R1_versor[2])

    df.plot(x='Distance', y='ForceI_x_expected', ax=ax, style='--')
    df.plot(x='Distance', y='ForceI_y_expected', ax=ax, style='-.')
    df.plot(x='Distance', y='ForceI_z_expected', ax=ax, style='--')
    plt.grid(visible=True)
    plt.xlabel('Distance R')
    plt.legend(
        [r'$Potential$', r'$F1_{x}$', r'$F1_{y}$', r'$F1_{z}$'])
    plt.title(handler.title+ r', $\nu = 1.6152500E-3$ - Expected')
    #plt.savefig('Figures/31_07/AxilrodTeller/Equilateral/ExpectedForce.png')
    plt.show()

def plot_computed_force(handler: GeometryHandler):
    df = pd.read_csv(file_location + handler.file_name, delimiter=',', index_col=False)

    ax = df.plot(x='Distance', y='Energy')
    df.plot(x='Distance', y='Force_x', ax=ax, style='--')
    df.plot(x='Distance', y='Force_y', ax=ax, style='-.')
    df.plot(x='Distance', y='Force_z', ax=ax, style='--')
    plt.grid(visible=True)
    plt.xlabel('Distance R')
    plt.legend(
        [r'$Potential$', r'$F1_{x}$', r'$F1_{y}$', r'$F1_{z}$'])
    plt.title(handler.title+ r', $\nu = 1.6152500E-3$ - Computed')
    #plt.savefig('Figures/31_07/AxilrodTeller/Equilateral/ComputedForce.png')
    plt.show()

def plot_potential(handler: GeometryHandler):
    df = pd.read_csv(file_location + handler.file_name, delimiter=',', index_col=False)

    ax = df.plot(x='Distance', y='Energy', color='blue', linestyle='-')

    nu = 1.6152500E-3
    def axilrodTellerEquilateral(x: float):
        return 11./ 8 * nu * x**-9
    x: np.array = np.linspace(0.3, 1.298999999999968, 1000)
    y: np.array = axilrodTellerEquilateral(x)
    ax.plot(x, y, color='orange', linestyle='--')

    plt.legend([r'$Potential$', r'$\frac{11}{8}\frac{ \nu }{R^9}$'])
    plt.xlabel('Distance R')
    plt.title(handler.title + r', $\nu = 1.6152500E-3$')
    #plt.savefig('Figures/31_07/AxilrodTeller/Equilateral/Potential.png')
    plt.show()

if __name__ == '__main__':
    handler = EquilateralGeometryHandler()
    plot_expected_force(handler)
    plot_computed_force(handler)
    plot_potential(handler)