import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.gridspec as gridspec
from matplotlib.lines import Line2D

file_location: str = '/Users/irene/TUM/thesis/AutoPas_Thesis/AutoPasFunctorBench/Argon/'


class Handler:
    def __init__(self):
        self.file_name: str = ''

        self.r1_versor: np.array = []

        self.title: str = ''
        self.xlabel: str = ''

        self.expected_force_fig: str = ''
        self.computed_force_fig: str = ''


class MovingAlongX(Handler):
    def __init__(self):
        super().__init__()
        self.file_name: str = file_location + 'Particle1AlongX.csv'

        self.r1_versor: np.array = [1, 0, 0]

        self.title: str = r'$r_1=[x,0,0]$, $r_2=[1,1,0]$, $r_3=[1,-1,0]$'
        self.xlabel: str = 'x'


class MovingAlongY(Handler):
    def __init__(self):
        super().__init__()
        self.file_name: str = file_location + 'Particle1AlongY.csv'

        self.r1_versor: np.array = [0, 1, 0]

        self.title: str = r'$r_1=[0,y,0]$, $r_2=[-1,1,0]$, $r_3=[1,1,0]$'
        self.xlabel: str = 'y'


class MovingAlongZ(Handler):
    def __init__(self):
        super().__init__()
        self.file_name: str = file_location + 'Particle1AlongZ.csv'

        self.r1_versor: np.array = [0, 0, 1]

        self.title: str = r'$r_1=[0,0,z]$, $r_2=[0,1,1]$, $r_3=[0,-1,1]$'
        self.xlabel: str = 'z'

color:dict = {'potential':'#006699', 'ForceX':'#669900', 'ForceY':'#990066', 'ForceZ':'#FF6600'}
def plot_subplot(handler: Handler, axs):
    df = pd.read_csv(handler.file_name, delimiter=',', index_col=False)

    #Plot computed
    df.plot(x='Position', y='Energy', ax=axs[0], color=color['potential'])
    df.plot(x='Position', y='Force_x', ax=axs[0], style='--', color=color['ForceX'])
    df.plot(x='Position', y='Force_y', ax=axs[0], style='-.', color=color['ForceY'])
    df.plot(x='Position', y='Force_z', ax=axs[0], style='--', color=color['ForceZ'])

    axs[0].set_xlabel(handler.xlabel)
    axs[0].set_title('Computed')
    axs[0].legend().set_visible(False)
    axs[0].grid(visible=True)

    #Plot expected
    df['dU_dr1'] = np.gradient(df['Energy'].values, df['Position'].values)
    df['ForceI_x_expected'] = -df['dU_dr1'] * handler.r1_versor[0]
    df['ForceI_y_expected'] = -df['dU_dr1'] * handler.r1_versor[1]
    df['ForceI_z_expected'] = -df['dU_dr1'] * handler.r1_versor[2]

    df.plot(x='Position', y='Energy', ax=axs[1], color=color['potential'])
    df.plot(x='Position', y='ForceI_x_expected', ax=axs[1], style='--', color=color['ForceX'])
    df.plot(x='Position', y='ForceI_y_expected', ax=axs[1], style='-.', color=color['ForceY'])
    df.plot(x='Position', y='ForceI_z_expected', ax=axs[1], style='--', color=color['ForceZ'])
    axs[1].set_xlabel(handler.xlabel)
    axs[1].set_title('Expected')
    axs[1].legend().set_visible(False)
    axs[1].grid(visible=True)

def plot():
    # gridspec inside gridspec
    fig = plt.figure(layout='constrained', figsize=(15, 10))
    subfigs = fig.subfigures(4, 1, wspace=0.07, height_ratios=[5,5,5,1])

    ## Moving Along X
    axsX = subfigs[0].subplots(1, 2, sharey=True)
    plot_subplot(MovingAlongX(), axsX)
    subfigs[0].suptitle(r'Moving Particle $i$ along $X$', fontsize='x-large')

    ## Moving Along Y
    axsY = subfigs[1].subplots(1, 2, sharey=True)
    plot_subplot(MovingAlongY(), axsY)
    subfigs[1].suptitle(r'Moving Particle $i$ along $Y$', fontsize='x-large')


    ## Moving Along Z
    axsZ = subfigs[2].subplots(1, 2, sharey=True)
    plot_subplot(MovingAlongZ(), axsZ)
    subfigs[2].suptitle(r'Moving Particle $i$ along $Z$', fontsize='x-large')

    fig.suptitle('Argon Potential', fontsize='xx-large')

    line_potential = Line2D([0], [0], label=r'$\Delta U$', color=color['potential'])
    line_forceX = Line2D([0], [0], label=r'$F_{i,x}$', color=color['ForceX'])
    line_forceY = Line2D([0], [0], label=r'$F_{i,y}$', color=color['ForceY'])
    line_forceZ = Line2D([0], [0], label=r'$F_{i,z}$', color=color['ForceZ'])
    handles = [line_potential, line_forceX, line_forceY, line_forceZ]
    axsLegend = subfigs[3].subplots(1,1)
    axsLegend.legend(handles=handles, ncol = 5, bbox_to_anchor=(0.4,1), loc='upper left')
    axsLegend.set_axis_off()

    plt.savefig('Figures/13_08/Argon/Complete.png')

    plt.show()


if __name__ == '__main__':
    plot()
