import numpy as np
import math
import matplotlib.pyplot as plt
from scipy.special import hermite
import tkinter as tk
from tkinter import ttk
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg



#constants
hbar = 1.0  
m = 1.0     
omega = 1.0 # Angular frequency
eV = 1.0    

initial_n = 0  #initial quantum number

N = 1000

def wavefunction_HO(n, x):
    prefactor = (m * omega / (np.pi * hbar)) ** 0.25
    norm = 1 / np.sqrt(2**n * math.factorial(n))
    # Hermite polynomial Hn
    hermite_poly = hermite(n)(np.sqrt(m * omega / hbar) * x)
    psi_n = prefactor * norm * hermite_poly * np.exp(- (m * omega * x**2) / (2 * hbar))
    return psi_n

def energy_HO(n):
    Ej = (n + 0.5) * hbar * omega
    EeV = Ej / eV
    return Ej, EeV

#create the main window
root = tk.Tk()
root.title("Quantum Harmonic Oscillator Wavefunction Explorer")
root.state('zoomed')

# Create figures with two subplots: Wavefunction and Probability Density
fig, (ax_wave, ax_prob) = plt.subplots(2, 1, figsize=(8, 6))

# Adjust layout for better spacing
plt.tight_layout(pad=3.0)

# Embed the figure into the Tkinter window
canvas = FigureCanvasTkAgg(fig, master=root)
canvas.draw()
canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=1)


def update_HO():
    # Get current quantum number from the slider
    n = int(round(slider_n.get()))
    
    # Define the x range (arbitrary units for visualization)
    x = np.linspace(-5, 5, N)
    
    # Calculate wavefunction and probability density
    psi_n = wavefunction_HO(n, x)
    prob_density = psi_n ** 2
    
    # Calculate energy levels
    Ej, EeV = energy_HO(n)
    
    # Clear previous plots
    ax_wave.cla()
    ax_prob.cla()
    
    # Plot wavefunction
    ax_wave.plot(x, psi_n, label=f"Wavefunction (n={n})", color='blue')
    ax_wave.set_title(f"Wavefunction for Quantum Harmonic Oscillator (n={n})")
    ax_wave.set_xlabel("Position x (arbitrary units)")
    ax_wave.set_ylabel(r'$\psi_n(x)$')
    ax_wave.axvline(x=0, color='black', linewidth=2, linestyle='--')  # Thicken x=0 axis
    ax_wave.axhline(y=0, color='black', linewidth=2, linestyle='--')  # Thicken y=0 axis
    ax_wave.grid(True)
    ax_wave.legend()
    
    # Plot probability density
    ax_prob.plot(x, prob_density, label=f"Probability Density (n={n})", color='red')
    ax_prob.set_title(f"Probability Density for Quantum Harmonic Oscillator (n={n})")
    ax_prob.set_xlabel("Position x (arbitrary units)")
    ax_prob.set_ylabel(r'$|\psi_n(x)|^2$')
    ax_prob.axvline(x=0, color='black', linewidth=2, linestyle='--')  # Thicken x=0 axis
    ax_prob.axhline(y=0, color='black', linewidth=2, linestyle='--')  # Thicken y=0 axis
    ax_prob.grid(True)
    ax_prob.legend()
    
    # Update energy label
    label_energy.config(text=f"Energy Level Eₙ = {EeV:.5f} eV")
    
    # Update current n label
    current_n.set(f"n = {n}")
    
    # Redraw the canvas to show updates
    canvas.draw()

# Frame to hold sliders and labels
frame_sliders = ttk.Frame(root)
frame_sliders.pack(side=tk.TOP, fill=tk.X, padx=20, pady=10)

# Label for quantum number slider
label_n = ttk.Label(frame_sliders, text="Quantum Number (n):")
label_n.pack(side=tk.LEFT, padx=(0, 10))

# Slider for quantum number n
slider_n = ttk.Scale(frame_sliders, from_=0, to=16, orient=tk.HORIZONTAL, command=lambda val: update_HO())
slider_n.set(initial_n)
slider_n.pack(side=tk.LEFT, padx=(0, 10))

# Display current n
current_n = tk.StringVar()
current_n.set(f"n = {initial_n}")
label_current_n = ttk.Label(frame_sliders, textvariable=current_n)
label_current_n.pack(side=tk.LEFT, padx=(0, 20))

# Frame to hold energy label
frame_energy = ttk.Frame(root)
frame_energy.pack(side=tk.TOP, fill=tk.X, padx=20, pady=10)

# Label to display energy
label_energy_text = ttk.Label(frame_energy, text="Energy Level E = ", font=("Arial", 14))
label_energy_text.pack(side=tk.LEFT)

label_energy = ttk.Label(frame_energy, text="", font=("Arial", 14, "bold"))
label_energy.pack(side=tk.LEFT)

# Function to update when slider is moved and released
def slidermove_HO(event):
    # Update n and call update_HO
    n = int(round(slider_n.get()))
    current_n.set(f"n = {n}")
    update_HO()

# Bind the slider to the callback function
slider_n.bind("<ButtonRelease-1>", slidermove_HO)



# Initial plot
update_HO()

# Start the Tkinter event loop
root.mainloop()
