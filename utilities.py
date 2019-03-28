import numpy as np

import matplotlib as mpl
import matplotlib.pyplot as plt

from ipywidgets import interact, fixed

mpl.rcParams["figure.figsize"] = (10, 6)
mpl.rcParams["legend.fontsize"] = 14

d21 = 7.37E-05 # in eV^2
d32 = 2.54E-03
d31 = d32 + d21

dm2 = [d21, d31, d32]

t12 = np.arcsin((1/3)**0.5)
t13 = np.arcsin(0.)
t23 = np.arcsin(0.5**0.5)

def PMNS(t12, t13, t23):
        s12 = np.sin(t12)
        c12 = np.cos(t12)
        s23 = np.sin(t23)
        c23 = np.cos(t23)
        s13 = np.sin(t13)
        c13 = np.cos(t13)
        return np.array([[ c12*c13, s12*c13, s13],
                                    [-s12*c23 - c12*s23*s13, c12*c23 - s12*s23*s13, s23],
                                    [ s12*s23 - c12*s23*s13,-c12*s23 - s12*c23*s13, c23]])



def plot_compare(prob_survival_exo4, posc, xlabel, x=None, ybounds=(0,1.)):
    def plot(L, E, t13=0, log=False, t12=t12, t23=t23, dm2=dm2):
        U = PMNS(t12, t13, t23)
        fig, ax = plt.subplots()
        
        if log:
            func = ax.semilogx
        else:
            func = ax.plot
            
        if x is None:
            x_ = LE
        else:
            x_ = x
        
        ax.set_ylim(0,1.35)
        ax.set_xlim(np.min(x), np.max(x))
        ax.set_xlabel(xlabel)
        ax.set_ylabel(r"Probability $\nu_e \rightarrow \nu_e$")
        func(x_, prob_survival_exo4(dm2[0], t12, L/E), label=f"E = {E} GeV\n"+r"$P_{ee} = 1 -\sin^2 2\theta_{12}\sin^2(\Delta m^2_{12}L/4E_\nu)$")
        func(x_, posc(0, 0, U, dm2, L/E), label="$P_{ee} = |<v_{e}|U|v_{e}>|^2$")
        plt.legend(loc="best")
        plt.grid()
        plt.show()
    return plot
    
def plot_posc(posc, xbounds, xlabel, x=None, ybounds=(0,1.)):
    def plot(L, E, t13=0, log=False, t12=t12, t23=t23, dm2=dm2):
        LE = L/E
        U = PMNS(t12, t13, t23)
        Pe = posc(0, 0, U, dm2, LE)
        Pm = posc(0, 1, U, dm2, LE)
        Pt = posc(0, 2, U, dm2, LE)
        fig, ax = plt.subplots()
                
        if x is None:
            x_ = LE
        else:
            x_ = x
                
        if log:
            func = ax.semilogx
        else:
            func = ax.plot
        
        func(x_, Pe, '-', label=r'$\nu_e \rightarrow \nu_e$')
        func(x_, Pm, 'k', label=r'$\nu_e\rightarrow\nu_\mu$')
        func(x_, Pt, 'r', label=r'$\nu_e\rightarrow\nu_\tau$')
        ax.set_xlabel(xlabel)
        ax.set_ylabel(r"Probability $\nu_e \rightarrow \nu_x$")
        ax.set_ylim(*ybounds)
        ax.set_xlim(*xbounds)
        plt.legend(loc="best")
        plt.grid()
        plt.show()
    return plot
    
def plot_posc_conv(posc, posc_conv, xbounds, xlabel, resolution, x=None, ybounds=(0,1.05), t13=0., dayabay=False):
    def plot(L, E, t13=t13, log=False, t12=t12, t23=t23, dm2=dm2):
        LE = L/E
        U0 = PMNS(t12, 0, t23)
        U = PMNS(t12, t13, t23)
        Pe = posc(0, 0, U, dm2, LE)
        Pe_conv = posc_conv(0, 0, U, dm2, LE, resolution)
        Pe_conv0 = posc_conv(0, 0, U0, dm2, LE, resolution)
        fig, ax = plt.subplots()
        
        if x is None:
            x_ = LE
        else:
            x_ = x
        
        if log:
            func = ax.semilogx
        else:
            func = ax.plot
            
        func(x_, Pe, '-', label=r'$\nu_e \rightarrow \nu_e$ without resolution', color="royalblue")
        func(x_, Pe_conv, 'k', label=r'$\nu_e\rightarrow\nu_e$ with resolution')
        func(x_, Pe_conv0, 'k--', label=r'$\nu_e\rightarrow\nu_e$ with resolution, $\theta_{13} = 0$')
    
        if dayabay:
            plt.axvline(x=0.47, linestyle="-.", color="tomato", label="near detector", linewidth=2.0)
            plt.axvline(x=1.65, linestyle="-.", color="orange", label="far detector", linewidth=2.0)
    
        ax.set_xlabel(xlabel)
        ax.set_ylabel(r"Probability $\nu_e \rightarrow \nu_x$")
        ax.set_ylim(*ybounds)
        ax.set_xlim(*xbounds)
        plt.legend(loc="best")
        plt.grid()
        plt.show()
    return plot
    
        
def interactive_plot(function, L=1.0, E=1.0, log=False):
    interact(function, L=fixed(L), E=fixed(E), 
             dm2=fixed(dm2), t13=(0.0,0.2, 0.001),
             t23=(0.0, 1.5, 0.001), t12=(0.0, 1.5, 0.001), log=log);
    