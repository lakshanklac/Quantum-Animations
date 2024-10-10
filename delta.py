import numpy as np
import matplotlib.pyplot as plt
import tkinter as tk
from tkinter import ttk
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg

hbar = 1  # reduced Planck's constant, J·s
me = 1    # mass of an electron, kg

N = 1000  # Number of points for plotting
initial_a = 0

# Define the delta function approximation (non-normalized)
def delta_approx(x, epsilon=0.01, a=1):
    # No normalization; just the raw exponential form
    delta = np.exp(-x**2 / (2 * epsilon**2)) 
    V = -a * delta  # Scale the delta function by 'a'
    return V

# Define the exponential decay function
def exponential_decay(x, a):
    B = np.sqrt(me * a)  # Scaling factor for exponential decay
    return B * np.exp(-np.abs(x))

# Create the main window
root = tk.Tk()
root.title("Quantum Harmonic Oscillator Wavefunction Explorer")
root.state('zoomed')

# Create a single figure for both plots
fig, plot_dirac = plt.subplots(figsize=(8, 6))

# Adjust layout for better spacing
plt.tight_layout(pad=3.0)

# Embed the figure into the Tkinter window
canvas = FigureCanvasTkAgg(fig, master=root)
canvas.draw()
canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=1)

def update_Dirac():
    # Get current quantum number from the slider
    a = int(round(slider_n.get()))

    # Define the x range (arbitrary units for visualization)
    x = np.linspace(-10, 10, N)

    # Calculate potential and wavefunction
    potential = delta_approx(x, a=a)  # Update the potential with 'a'
    psi_n = exponential_decay(x, a)  # Update the wavefunction with 'a'
    mag = -potential

    # Clear previous plots
    plot_dirac.cla()

    # Plot both potential and wavefunction in the same graph
    plot_dirac.plot(x, mag, linestyle=':', color='b', lw=2, label='Magnitude of Potential')
    plot_dirac.plot(x, potential, color='r', lw=2, label=r'$V(x) = -a*\delta(x)$')
    plot_dirac.plot(x, psi_n, color='g', lw=2, label=r'$y = B e^{-|x|}$')

    # Set titles and labels
    plot_dirac.set_title('Delta Function Potential and Exponential Decay')
    plot_dirac.set_xlabel('Position x (arbitrary units)')
    plot_dirac.set_ylabel('y, V(x)')
    plot_dirac.legend()
    plot_dirac.grid(True)

    # Update current 'a' label
    current_n.set(f"a = {a}")

    # Redraw the canvas to show updates
    canvas.draw()

# Frame to hold sliders and labels
frame_sliders = ttk.Frame(root)
frame_sliders.pack(side=tk.TOP, fill=tk.X, padx=20, pady=10)

# Label for constant a slider
label_n = ttk.Label(frame_sliders, text="Constant (a):")
label_n.pack(side=tk.LEFT, padx=(0, 10))

# Slider for constant a
slider_n = ttk.Scale(frame_sliders, from_=0, to=100, orient=tk.HORIZONTAL, command=lambda val: update_Dirac())
slider_n.set(initial_a)
slider_n.pack(side=tk.LEFT, padx=(0, 10))

# Display current a
current_n = tk.StringVar()
current_n.set(f"a = {initial_a}")
label_current_n = ttk.Label(frame_sliders, textvariable=current_n)
label_current_n.pack(side=tk.LEFT, padx=(0, 20))

# Function to update when slider is moved and released
def slidermove_dirac(event):
    # Update a and call update_Dirac
    a = float(round(slider_n.get()))
    current_n.set(f"a = {a}")
    update_Dirac()

# Bind the slider to the callback function
slider_n.bind("<ButtonRelease-1>", slidermove_dirac)

# Initial plot
update_Dirac()

# Start the Tkinter event loop
root.mainloop()
