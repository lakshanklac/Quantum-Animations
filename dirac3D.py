import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import tkinter as tk
from tkinter import ttk
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg

hbar = 1  # reduced Planck's constant
me = 1    # mass of an electron

N = 100  # Number of points for plotting (reduced for faster plotting)
initial_a = 1

# Define the delta function approximation (non-normalized)
def delta_approx(x, y, epsilon=0.5, a=1):
    delta = np.exp(-(x**2 + y**2) / (2 * epsilon**2))
    V = -a * delta
    return V

# Define the exponential decay function
def exponential_decay(x, y, a):
    B = np.sqrt(me * a)
    r = np.sqrt(x**2 + y**2)  # Radial distance
    return B * np.exp(-r)

# Create the main window
root = tk.Tk()
root.title("3D Delta Function Well")
root.state('zoomed')

# Create a figure for the 3D plot
fig = plt.figure(figsize=(8, 6))
plot_dirac = fig.add_subplot(111, projection='3d')

# Embed the figure into the Tkinter window
canvas = FigureCanvasTkAgg(fig, master=root)
canvas.draw()
canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=1)

# Function to update the 3D plot based on slider value
def update_Dirac():
    a = slider_n.get()

    # Create a grid in both x and y directions
    x = np.linspace(-10, 10, N)
    y = np.linspace(-10, 10, N)
    X, Y = np.meshgrid(x, y)

    # Calculate the potential and wavefunction
    potential = delta_approx(X, Y, a=a)
    psi_n = exponential_decay(X, Y, a)
    
    # Clear previous plots
    plot_dirac.cla()

    # Plot the potential (Delta Function Well)
    plot_dirac.plot_surface(X, Y, potential, alpha=0.7, cmap='Reds', rstride=1, cstride=1)
    
    # Plot the wavefunction (Exponential Decay)
    plot_dirac.plot_surface(X, Y, psi_n, alpha=0.6, cmap='Blues', rstride=1, cstride=1)
    
    # Set titles and labels
    plot_dirac.set_title('3D Delta Function Well and Exponential Decay')
    plot_dirac.set_xlabel('X axis')
    plot_dirac.set_ylabel('Y axis')
    plot_dirac.set_zlabel('Potential / Wavefunction')
    
    # Update the canvas with the new plot
    canvas.draw()

# Frame to hold sliders and labels
frame_sliders = ttk.Frame(root)
frame_sliders.pack(side=tk.TOP, fill=tk.X, padx=20, pady=10)

# Label for constant a slider
label_n = ttk.Label(frame_sliders, text="Constant (a):")
label_n.pack(side=tk.LEFT, padx=(0, 10))

# Slider for constant a
slider_n = ttk.Scale(frame_sliders, from_=0.1, to=5, orient=tk.HORIZONTAL, length=400, command=lambda val: update_Dirac())
slider_n.set(initial_a)
slider_n.pack(side=tk.LEFT, padx=(0, 10))

# Display current value of a
current_n = tk.StringVar()
current_n.set(f"a = {initial_a}")
label_current_n = ttk.Label(frame_sliders, textvariable=current_n)
label_current_n.pack(side=tk.LEFT, padx=(0, 20))

# Function to update when slider is moved
def slidermove_dirac(event):
    a = float(slider_n.get())
    current_n.set(f"a = {a:.2f}")
    update_Dirac()

# Bind the slider to update on release
slider_n.bind("<ButtonRelease-1>", slidermove_dirac)

# Initial plot
update_Dirac()

# Start the Tkinter event loop
root.mainloop()
