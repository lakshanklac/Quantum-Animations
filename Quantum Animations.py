import numpy as np
import matplotlib.pyplot as plt
import tkinter as tk
from tkinter import ttk
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg

#constants
hbar=1.054571817e-34 #reduced Planck's constant,J·s
me=9.10938356e-31 #mass of an electron,kg
eV=1.602176634e-19 #electron volt,J

#initial parameters
initial_n=1 #quantum number
initial_L_nm=1.0 #initial well width in nm

#Convert initial_L to meters
initial_L=initial_L_nm * 1e-9  #meters

#Spatial grid parameters
N=10000  #Number of spatial points


#function to calculate wavefunction_infinitesquarewell
def wavefunction_infinitesquarewell(n, x, L):
    y=np.sqrt(2 / L) * np.sin(n * np.pi * x / L)
    return y

#function to calculate energy_infinitesquarewell
def energy_infinitesquarewell(n, L):
    EJ=(n**2 * (np.pi**2) * (hbar**2)) / (2 * me * (L**2))
    EeV=EJ / eV
    return EJ, EeV


#create the main window
root=tk.Tk()
root.title("Infinite Square Well wavefunction_infinitesquarewell Explorer")
root.state('zoomed')
#create figures
fig, (plotwave, plotdensity)=plt.subplots(2, 1, figsize=(8, 6))

#gap between two plots
plt.tight_layout(pad=3.0)

canvas=FigureCanvasTkAgg(fig, master=root)
canvas.draw()
canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=1)

#function to update_infinitesquarewell through the sliders
def update_infinitesquarewell():
    #get current values from sliders
    n=int(round(slider_n.get()))
    L_nm=float(slider_L.get())
    L=L_nm * 1e-9

 
    x=np.linspace(0, L, N)  #meters

    psi_n=wavefunction_infinitesquarewell(n, x, L)
    prob_density=psi_n**2

    #re-calculate energy_infinitesquarewell
    EJ, EeV=energy_infinitesquarewell(n, L)

    #clear previous plots
    plotwave.cla()
    plotdensity.cla()

    #plot wavefunction_infinitesquarewell
    plotwave.plot(x * 1e9, psi_n, color='blue', label=f'n={n}')
    plotwave.set_title("wavefunction_infinitesquarewell in Infinite Square Well")
    plotwave.set_xlabel("Position x (nm)")
    plotwave.set_ylabel(r'$\psi_n(x)$')
    plotwave.grid(True)
    plotwave.axvline(x=0, color='black', linestyle='--')
    plotwave.axvline(x=L * 1e9, color='black', linestyle='--')
    plotwave.legend()
    
    #thicken the x=0 axis and y=0 axis for the wavefunction_infinitesquarewell plot
    plotwave.axvline(x=0, color='black', linestyle='--', linewidth=2)  # x=0
    plotwave.axhline(y=0, color='black', linewidth=2)  # y=0
    plotwave.axvline(x=L * 1e9, color='black', linestyle='--')

    #plot probability density
    plotdensity.plot(x * 1e9, prob_density, color='red', label=f'n={n}')
    plotdensity.set_title("Probability Density in Infinite Square Well")
    plotdensity.set_xlabel("Position x (nm)")
    plotdensity.set_ylabel(r'$|\psi_n(x)|^2$')
    plotdensity.grid(True)
    plotdensity.axvline(x=0, color='black', linestyle='--')
    plotdensity.axvline(x=L * 1e9, color='black', linestyle='--')
    plotdensity.legend()
    
    #thicken the x=0 axis and y=0 axis for the probability density plot
    plotdensity.axvline(x=0, color='black', linestyle='--', linewidth=2)  # x=0
    plotdensity.axhline(y=0, color='black', linewidth=2)  # y=0
    plotdensity.axvline(x=L * 1e9, color='black', linestyle='--')


    #update_infinitesquarewell energy_infinitesquarewell label
    label_energy_infinitesquarewell.config(text=f"energy_infinitesquarewell Level E={EeV:.5f} eV")

    #update_infinitesquarewell current n and L labels
    current_n.set(f"n={n}")
    current_L.set(f"L={L_nm:.1f} nm")

    #r the canvas
    canvas.draw()


#frame to hold sliders
frame_sliders=ttk.Frame(root)
frame_sliders.pack(side=tk.TOP, fill=tk.X, padx=20, pady=10)

#label for quantum number slider
label_n=ttk.Label(frame_sliders, text="Quantum Number (n):")
label_n.pack(side=tk.LEFT, padx=(0,10))

#slider for quantum number n
slider_n=ttk.Scale(frame_sliders, from_=1, to=10, orient=tk.HORIZONTAL, command=lambda val: update_infinitesquarewell())
slider_n.set(initial_n)
slider_n.pack(side=tk.LEFT, padx=(0,10))

#display current n
current_n=tk.StringVar()
current_n.set(f"n={initial_n}")
label_current_n=ttk.Label(frame_sliders, textvariable=current_n)
label_current_n.pack(side=tk.LEFT, padx=(0,20))


#label for well width slider
label_L=ttk.Label(frame_sliders, text="Well Width (L) [nm]:")
label_L.pack(side=tk.LEFT, padx=(0,10))

#slider for well width L
slider_L=ttk.Scale(frame_sliders, from_=0.1, to=10.0, orient=tk.HORIZONTAL, command=lambda val: update_infinitesquarewell())
slider_L.set(initial_L_nm)
slider_L.pack(side=tk.LEFT, padx=(0,5))

#display current L
current_L=tk.StringVar()
current_L.set(f"L={initial_L_nm:.1f} nm")
label_current_L=ttk.Label(frame_sliders, textvariable=current_L)
label_current_L.pack(side=tk.LEFT, padx=(0,20))

#frame to hold energy_infinitesquarewell label
frame_energy_infinitesquarewell=ttk.Frame(root)
frame_energy_infinitesquarewell.pack(side=tk.TOP, fill=tk.X, padx=20, pady=10)

#dynamic label to display energy_infinitesquarewell
label_energy_infinitesquarewell=ttk.Label(frame_energy_infinitesquarewell, text="", font=("Arial", 14, "bold"))
label_energy_infinitesquarewell.pack(side=tk.LEFT)

#function to update_infinitesquarewell when sliders are moving
def slidermove_infinitesquarewell(event):
    #update_infinitesquarewell current_n and current_L labels
    n=int(round(slider_n.get()))
    current_n.set(f"n={n}")
    L_nm=float(slider_L.get())
    current_L.set(f"L={L_nm:.1f} nm")
    #update_infinitesquarewell plots
    update_infinitesquarewell()

#bind the sliders to the callback function
slider_n.bind("<ButtonRelease-1>", slidermove_infinitesquarewell)
slider_L.bind("<ButtonRelease-1>", slidermove_infinitesquarewell)

update_infinitesquarewell()

root.mainloop()
