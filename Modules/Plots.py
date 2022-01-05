from matplotlib import pyplot as plt
import numpy as np
import pandas as pd


def plot_bound(m_a, g_a, x_range, y_range, file_name, label_name, color_plot):
    """
        Plots the results of the analysis in a figure like Fig. 5 of arXiv:2201.XXXX .

        Parameters
        ----------
        m_a : array_like
            The masses of the axion which the limits apply to, in units of 10^(-22) eV.
        g_a : array_like
            The lower limit on the coupling g_a\gamma, in units of GeV^(-1).
        x_range : array_like
            The range of the x-axis (m_a, in units of 10^(-22) eV) that is to be displayed in the figure.
        y_range : array_like
            The range of the y-axis (g_a, in units of GeV^(-1)) that is to be displayed in the figure.
        file_name : string
            The name of the file to save the figure.
        label_name : string
            The name of the line in the figure corresponding to the new result to be plotted.
        color_plot : string
            The color of the line in the figure corresponding to the new result to be plotted.

        Returns
        -------
        None

        Notes
        -----
        This function is used to plot the results of the lower limits on the Ultra-Light axion coupling to the photon.

    """
    df_AGN_VLBA = pd.read_csv('InputData/PastBounds/AGN_VLBA.dat',
                      delim_whitespace=True, header=None)
    df_Keck_Low = pd.read_csv(
    'InputData/PastBounds/KeckLow.dat', delim_whitespace=True, header=None)
    df_Keck_High = pd.read_csv(
    'InputData/PastBounds/KeckHigh.dat', delim_whitespace=True, header=None)
    df_Castillo = pd.read_csv('InputData/PastBounds/Castilloetal.txt',delim_whitespace=True, header=None)

    g_Castillo  = df_Castillo[1]
    m_Castillo  = df_Castillo[0]

    g_AGN_VLBA = df_AGN_VLBA[1]
    m_AGN_VLBA = df_AGN_VLBA[0]

    g_Keck_Low = df_Keck_Low[1]
    m_Keck_Low = df_Keck_Low[0]
    m_Keck_Low = [10**22*i for i in m_Keck_Low]

    g_Keck_High = df_Keck_High[1]
    m_Keck_High = df_Keck_High[0]
    m_Keck_High = [10**22*i for i in m_Keck_High]

    fig, ax = plt.subplots()

    ### Your results ###
    plt.plot(m_a, g_a, '-', color=color_plot)
    plt.fill_between(m_a, g_a, 1, color=color_plot, alpha=.3, label=label_name)

    ##### Bounds to compare ######
    ### AGNS arXiv:1811.10997###
    plt.plot(m_AGN_VLBA, g_AGN_VLBA, '--', alpha=.8, color='orangered')
    plt.fill_between(m_AGN_VLBA, g_AGN_VLBA, 1, color='orangered', alpha=.1, label='VLBA')

    ## CAST arXiv:1705.02290##
    plt.plot([.01, 1000], [0.66*pow(10, -10), 0.66*pow(10, -10)],
            ls='-.', color='darkgrey')  # CAST
    plt.annotate(r'CAST', xy=(450, 0.00000000007), fontsize=22)

    ## SN1987 arXiv:1410.3747##

    plt.plot([.01, 1000], [5.3*pow(10, -12), 5.3 * 
            pow(10, -12)],  color='darkgrey')  # sn1987
    plt.annotate(r'SN $1987$A', xy=(325, 0.0000000000056), fontsize=22)
    
    ## Keck-BICEP arXiv:2108.03316##

    plt.plot(m_Keck_Low, g_Keck_Low, '--', alpha=.8, color='forestgreen')
    plt.fill_between(m_Keck_Low, g_Keck_Low, 1)
    plt.plot(m_Keck_High, g_Keck_High, '--', alpha=.8, color='forestgreen')
    plt.fill_between(m_Keck_High, g_Keck_High, 1,
                    color='forestgreen', alpha=.1, label='Keck-BICEP')


    ## A. Castillo et al. ##
    plt.plot(m_Castillo, g_Castillo, '--', alpha= .8, color='mediumblue')
    plt.fill_between(m_Castillo, g_Castillo, 1, color='royalblue',
                    alpha=.3, label='A. Castillo et al.')



    ### Axis labels and range ###
    ax.set_xlabel(r'$m_a$ [$10^{-22}$ eV]', fontsize=34)
    ax.set_ylabel(r'$g_{a \gamma}$ [GeV$^{-1}]$', fontsize=34)

    ax.set_ylim([y_range[0], y_range[1]])
    ax.set_xlim([x_range[0], x_range[1]])

    ### Ticks parameters ####
    plt.tick_params(labelbottom=True, labeltop=False, labelleft=True, labelright=False,
                bottom=True, top=True, left=True, right=True, direction='in', which='major', pad=5.5, width=1.2, length=13, labelsize=26)

    plt.tick_params(labelbottom=True, labeltop=False, labelleft=True, labelright=False,
                bottom=True, top=True, left=True, right=True, direction='in', which='minor', pad=5.5, width=1.2, length=10, labelsize=26)

    ### Plot legend ###
    plt.legend()

    ### Save the figure ###
    plt.savefig('Output/Figures/'+file_name+'.pdf')

