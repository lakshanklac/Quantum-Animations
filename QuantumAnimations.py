import tkinter as tk
from tkinter import ttk
import numpy as np
import matplotlib.pyplot as plt
from scipy.special import hermite
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import math
from scipy.special import lpmv
from mpl_toolkits.mplot3d import Axes3D
from tkinter import messagebox
from scipy.special import legendre
from scipy.special import sph_harm


# Function to open a new window for Button Infinite Square Well
def open_window_infinitesquarewell():
    #constants
    hbar=1.054571817e-34 #reduced Planck's constant,J·s
    me=9.10938356e-31 #mass of an electron,kg
    eV=1.602176634e-19 #electron volt,J
    
    #initial parameters
    initial_n=1 #quantum number
    initial_L_nm=0.1 #initial well width in nm
    
    #Convert initial_L to meters
    #initial_L=initial_L_nm * 1e-9  #meters
    
    N=1000  #Number of spatial points
    
    
    # Function to calculate wavefunction_infinitesquarewell
    def wavefunction_infinitesquarewell(n, x, L):
        y=np.sqrt(2 / L) * np.sin(n * np.pi * x / L)
        return y
    
    # Function to calculate energy_infinitesquarewell
    def energy_infinitesquarewell(n, L):
        EJ=(n**2 * (np.pi**2) * (hbar**2)) / (2 * me * (L**2))
        EeV=EJ / eV
        return EJ, EeV
    
    
    # Create the main window
    root=tk.Tk()
    root.title("Infinite Square Well")
    root.state('zoomed')
    
    # Create figures
    fig, (plotwave, plotdensity)=plt.subplots(2, 1, figsize=(8, 6))
    
    # Gap between two plots
    plt.tight_layout(pad=3.0)
    
    canvas=FigureCanvasTkAgg(fig, master=root)
    canvas.draw()
    canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=1)
    
    # Function to update_infinitesquarewell through the sliders
    def update_infinitesquarewell():
        #get current values from sliders
        n=int(round(slider_n.get()))
        L_nm=float(slider_L.get())
        L=L_nm * 1e-9
    
     
        x=np.linspace(0, L, N)
    
        psi_n=wavefunction_infinitesquarewell(n, x, L)
        prob_density=psi_n**2
    
        # Re-calculate energy_infinitesquarewell
        EJ, EeV=energy_infinitesquarewell(n, L)
    
        #clear previous plots
        plotwave.cla()
        plotdensity.cla()
    
        L=float(L)
        #plot wavefunction_infinitesquarewell
        plotwave.plot(x * 1e9, psi_n, color='blue', label=f'n={n}')
        plotwave.set_title(f"Wavefunction of Infinite Square Well for (n={n}) and (L={L:.2e})")
        plotwave.set_xlabel("x (nm)")
        plotwave.set_ylabel(r'$\psi_n(x)$')
        plotwave.grid(True)
        plotwave.axvline(x=0, color='black', linestyle='--')
        plotwave.axvline(x=L * 1e9, color='black', linestyle='--')
        plotwave.legend()
        
        #thicken the x=0 axis and y=0 axis for the wavefunction_infinitesquarewell plot
        plotwave.axvline(x=0, color='black', linestyle='--')  # x=0
        plotwave.axhline(y=0, color='black', linestyle='--')  # y=0
        plotwave.axvline(x=L * 1e9, color='black', linestyle='--')
    
        #plot probability density
        plotdensity.plot(x * 1e9, prob_density, color='red', label=f'n={n}')
        plotdensity.set_title(f'Probability Density of Infinite Square Well for (n={n}) and (L={float(L):.2e})')
        plotdensity.set_xlabel("x (nm)")
        plotdensity.set_ylabel(r'$|\psi_n(x)|^2$')
        plotdensity.grid(True)
        plotdensity.axvline(x=0, color='black', linestyle='--')
        plotdensity.axvline(x=L * 1e9, color='black', linestyle='--')
        plotdensity.legend()
        
        #thicken the x=0 axis and y=0 axis for the probability density plot
        plotdensity.axvline(x=0, color='black', linestyle='--')  # x=0
        plotdensity.axhline(y=0, color='black', linestyle='--')  # y=0
        plotdensity.axvline(x=L * 1e9, color='black', linestyle='--')
    
    
        #update_infinitesquarewell energy_infinitesquarewell label
        label_energy_infinitesquarewell.config(text=f"Energy of Level n={n:.0f} Level E={EeV:.3f} eV")
    
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
    slider_n=ttk.Scale(frame_sliders, from_=1, to=10, orient=tk.HORIZONTAL,length=500, command=lambda val: update_infinitesquarewell())
    slider_n.set(initial_n)
    slider_n.pack(side=tk.LEFT, padx=(0,10))
    
    #display current n
    current_n=tk.StringVar()
    current_n.set(f"n={initial_n}")
    label_current_n=ttk.Label(frame_sliders, textvariable=current_n)
    label_current_n.pack(side=tk.LEFT, padx=(0,20))
    
    
    #label for well width slider
    label_L=ttk.Label(frame_sliders, text="Well Width (L) (nm):")
    label_L.pack(side=tk.LEFT, padx=(0,10))
    
    #slider for well width L
    slider_L=ttk.Scale(frame_sliders, from_=0.1, to=10.0, orient=tk.HORIZONTAL,length=500, command=lambda val: update_infinitesquarewell())
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
    
 

# Function to open a new window for Button B
def open_window_HO():
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
        ax_wave.set_xlabel("x (arbitrary units)")
        ax_wave.set_ylabel(r'$\psi_n(x)$')
        ax_wave.axvline(x=-4, color='r',linestyle='--')
        ax_wave.axvline(x=-2, color='black',linestyle='--')
        ax_wave.axhline(y=0, color='black',linestyle='--')# Thicken x=0 axis
        ax_wave.axvline(x=2, color='black', linestyle='--')  # Thicken y=0 axis
        ax_wave.axvline(x=4, color='r',linestyle='--')
        ax_wave.grid(True)
        ax_wave.legend()
        
 
        
        # Plot probability density
        ax_prob.plot(x, prob_density, label=f"Probability Density (n={n})", color='red')
        ax_prob.set_title(f"Probability Density for Quantum Harmonic Oscillator (n={n})")
        ax_prob.set_xlabel("Position x (arbitrary units)")
        ax_prob.set_ylabel(r'$|\psi_n(x)|^2$')
        ax_prob.axvline(x=-4, color='b',  linestyle='--')
        ax_prob.axvline(x=-2, color='black',  linestyle='--')  # Thicken x=0 axis
        ax_prob.axhline(y=0, color='black',  linestyle='--')
        ax_prob.axvline(x=2, color='black', linestyle='--')
        ax_prob.axvline(x=4, color='b', linestyle='--')# Thicken y=0 axis
        ax_prob.grid(True)
        ax_prob.legend()
        
        energy_factor=(n+0.5)
        # Update energy label
        hbar_symbol = '\u210F'  # Unicode for ħ
        omega_symbol = '\u03C9'  # Unicode for ω (omega)
        label_energy.config(text=f"Energy Level E_{n:.0f} = {energy_factor:.2f} * {hbar_symbol} * {omega_symbol}")


        
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
    slider_n = ttk.Scale(frame_sliders, from_=0, to=16, orient=tk.HORIZONTAL,length=500, command=lambda val: update_HO())
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
    #label_energy_text = ttk.Label(frame_energy, text="Energy Level E = ", font=("Arial", 14))
    #label_energy_text.pack(side=tk.LEFT)
    
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
    
def open_window_delta_dirac_2D():
    
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
        
        alpha_symbol = '\u03B1'  # Unicode for alpha

        # Plot both potential and wavefunction in the same graph
        alpha_symbol = '\u03B1'  # Unicode for alpha

        # Plot both potential and wavefunction in the same graph
        plot_dirac.plot(x, mag, linestyle=':', color='b', lw=2, label='Magnitude of Potential')
        plot_dirac.plot(x, potential, color='r', lw=2, label=fr'$V(x) = -{alpha_symbol}*\delta(x)$')
        plot_dirac.plot(x, psi_n, color='g', lw=2, label=fr'$y = \sqrt{{m{alpha_symbol}}}/{{\hbar}} {alpha_symbol} e^{{-\left(\frac{{m{alpha_symbol}}}{{\hbar^2}}\right)|x|}}$')
                                    
    
    
        # Set titles and labels
        plot_dirac.set_title('Bound State Wave Function of Delta Function Potential Well')
        plot_dirac.set_xlabel('x (arbitrary units)')
        plot_dirac.set_ylabel(r"$\psi(x)$ , $V(x)$")
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
    slider_n = ttk.Scale(frame_sliders, from_=0, to=100, orient=tk.HORIZONTAL,length=500, command=lambda val:     update_Dirac())    
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
    
def open_window_delta_scatter():
    hbar = 1  # reduced Planck's constant, J·s
    me = 1    # mass of an electron, kg
    
    N = 1000  # Number of points for plotting
    initial_a = 10
    initial_E = 10
    
    # Define the delta function approximation (non-normalized)
    def delta_approx(x, epsilon=0.01, a=1):
        delta = np.exp(-x**2 / (2 * epsilon**2))  # No normalization; just the raw exponential form
        V = -a * delta  # Scale the delta function by 'a'
        return V
    
    # Calculate the wavefunction for scattering states
    def calculate_wavefunctions(a, E):
        k = np.sqrt(2 * E)  # wavevector from energy
        beta = a / np.sqrt(2 * E)  # scattering strength related to the delta potential
    
        # Define regions for x < 0 and x > 0
        x_neg = np.linspace(-10, 0, N//2)
        x_pos = np.linspace(0, 10, N//2)
    
        A = 1
        B = 1j * beta * A / (1 - 1j * beta)
        F = A / (1 - 1j * beta)
    
        # Compute the wavefunctions for x < 0 and x > 0
        psi_A = A * np.exp(1j * k * x_neg)
        psi_B = B * np.exp(-1j * k * x_neg)
        psi_F = F * np.exp(1j * k * x_pos)
    
        return x_neg, x_pos, psi_A, psi_B, psi_F
    
    # Create the main window
    root = tk.Tk()
    root.title("Scattering State of Dirac Delta Potential")
    root.state('zoomed')
    
    # Create a single figure for both plots
    fig, plot_dirac = plt.subplots(figsize=(8, 6))
    
    # Adjust layout for better spacing
    plt.tight_layout(pad=3.0)
    
    def update_Dirac():
        # Get current 'a' and 'E' from the sliders
        a = float(slider_a.get())
        E = float(slider_E.get())
    
        # Define the x range and potential
        x = np.linspace(-10, 10, N)
        potential = delta_approx(x, a=a)  # Update the potential with 'a'
        mag = -potential
    
        # Get the updated wavefunctions based on 'a' and 'E'
        x_neg, x_pos, psi_A, psi_B, psi_F = calculate_wavefunctions(a, E)
    
        # Clear previous plots
        plot_dirac.cla()
    
        # Calculate Reflection (R) and Transmission (T) coefficients
        gamma = 2 * hbar**2 * E / (me * a**2)
        R = 1 / (1 + gamma)  # Reflection coefficient
        T = 1 - R            # Transmission coefficient (for conservation of probability)
    
        # Plot the delta potential
        plot_dirac.plot(x, mag, linestyle=':', color='b', lw=2, label='Magnitude of Potential')
        plot_dirac.plot(x, potential, color='r', lw=2, label=r'$V(x) = -a*\delta(x)$')
    
        # Plot the reflection and transmission as constant lines
        plot_dirac.axhline(y=R, linestyle='--', color='cyan', label='R (Reflection Coefficient)')
        plot_dirac.axhline(y=T, linestyle='--', color='darkblue', label='T (Transmission Coefficient)')
    
        # Add text annotations for R and T
        plot_dirac.text(0, R, f'R = {R:.4f}', fontsize=10, color='cyan', ha='left')
        plot_dirac.text(0, T, f'T = {T:.4f}', fontsize=10, color='darkblue', ha='left')
    
        # Plot the scattering wavefunctions
        plot_dirac.plot(x_neg, np.imag(psi_A), label=r"Im($\psi_A(x)$) for $x < 0$", color="b")
        plot_dirac.plot(x_neg, np.imag(psi_B), label=r"Im($\psi_B(x)$) for $x < 0$", color="r", linestyle="--")
        plot_dirac.plot(x_pos, np.imag(psi_F), label=r"Im($\psi_F(x)$) for $x > 0$", color="orange", linestyle="--")
    
        # Set titles and labels
        plot_dirac.set_title('Scattering State of Dirac Delta Potential')
        plot_dirac.set_xlabel('Position x (arbitrary units)')
        plot_dirac.set_ylabel(r"$\psi(x)$ and $V(x)$")
        plot_dirac.legend()
        plot_dirac.grid(True)
    
        # Update current 'a' and 'E' labels
        current_a.set(f"a = {a:.2f}")
        current_E.set(f"E = {E:.2f}")
    
        # Redraw the canvas to show updates
        canvas.draw()
    
    # Embed the figure into the Tkinter window
    canvas = FigureCanvasTkAgg(fig, master=root)
    canvas.draw()
    canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=1)
    
    # Frame to hold sliders and labels
    frame_sliders = ttk.Frame(root)
    frame_sliders.pack(side=tk.TOP, fill=tk.X, padx=20, pady=10)
    
    # Label for constant 'a' slider
    label_a = ttk.Label(frame_sliders, text="Potential strength (a):")
    label_a.pack(side=tk.LEFT, padx=(0, 10))
    
    # Slider for constant 'a'
    slider_a = ttk.Scale(frame_sliders, from_=0.1, to=20.0, orient=tk.HORIZONTAL,length=500, command=lambda val: update_Dirac())
    slider_a.set(initial_a)
    slider_a.pack(side=tk.LEFT, padx=(0, 10))
    
    # Display current 'a'
    current_a = tk.StringVar()
    current_a.set(f"a = {initial_a:.2f}")
    label_current_a = ttk.Label(frame_sliders, textvariable=current_a)
    label_current_a.pack(side=tk.LEFT, padx=(0, 20))
    
    # Label for energy 'E' slider
    label_E = ttk.Label(frame_sliders, text="Energy (E):")
    label_E.pack(side=tk.LEFT, padx=(0, 10))
    
    # Slider for energy 'E'
    slider_E = ttk.Scale(frame_sliders, from_=0.1, to=20.0, orient=tk.HORIZONTAL,length=500, command=lambda val: update_Dirac())
    slider_E.set(initial_E)
    slider_E.pack(side=tk.LEFT, padx=(0, 10))
    
    # Display current 'E'
    current_E = tk.StringVar()
    current_E.set(f"E = {initial_E:.2f}")
    label_current_E = ttk.Label(frame_sliders, textvariable=current_E)
    label_current_E.pack(side=tk.LEFT, padx=(0, 20))
    
    # Function to update when slider is moved and released
    def slidermove_dirac(event):
        update_Dirac()
    
    # Bind the sliders to the callback function
    slider_a.bind("<ButtonRelease-1>", slidermove_dirac)
    slider_E.bind("<ButtonRelease-1>", slidermove_dirac)
    
    # Initial plot
    update_Dirac()
    
    # Start the Tkinter event loop
    root.mainloop()

def open_window_legendre_poly():    
    
    # Create the main window
    root = tk.Tk()
    root.title("Legendre Polynomials")
    root.state('zoomed')
    
    # Create figures with two subplots
    fig, (ply) = plt.subplots(figsize=(8, 6))
    
    # Gap between two plots
    plt.tight_layout(pad=3.0)
    
    canvas = FigureCanvasTkAgg(fig, master=root)
    canvas.draw()
    canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=1)
    
    initial_n = 0
    
    # Function to update the Legendre polynomial plot through the sliders
    def update_legendre():
        # Get current values from sliders
        n = int(round(slider_n.get()))
        x = np.linspace(-1, 1, 4000)
    
        # Clear previous plots
        ply.cla()        # Clear the first subplot for Legendre polynomial
    
        P_n = legendre(n)
    
        # Evaluate P_n on the x values
        y = P_n(x)
    
        # Plot the polynomial on the specified axes
        ply.plot(x, y, label=f'(l={n})', lw=2)
    
    
        # Setting labels and title for the first plot
        ply.set_xlabel('x')
        ply.set_ylabel('P_l')
        ply.set_title(f"Legendre Polynomials of (l={n})" )
        ply.axhline(0, color='black', lw=0.5, ls='--')
        ply.axhline(y=0, color='black',ls='--')
        ply.axvline(0, color='black', lw=0.5, ls='--')
        ply.grid(True)
        ply.legend()
    
    
        
        # Update current n label
        current_n.set(f"n={n}")
    
        # Redraw the canvas to show updates
        canvas.draw()
    
    # Frame to hold sliders
    frame_sliders = ttk.Frame(root)
    frame_sliders.pack(side=tk.TOP, fill=tk.X, padx=20, pady=10)
    
    # Label for quantum number slider
    label_n = ttk.Label(frame_sliders, text="Azimuthal Quantum Number (l):")
    label_n.pack(side=tk.LEFT, padx=(0, 10))
    
    # Slider for quantum number n
    slider_n = ttk.Scale(frame_sliders, from_=0, to=50, orient=tk.HORIZONTAL,length=500, command=lambda val: update_legendre())
    slider_n.set(initial_n)
    slider_n.pack(side=tk.LEFT, padx=(0, 10))
    
    # Display current n
    current_n = tk.StringVar()
    label_current_n = ttk.Label(frame_sliders, textvariable=current_n)
    label_current_n.pack(side=tk.LEFT, padx=(0, 20))
    
    # Function to update when sliders are moving
    def slidermove_legendrepoly(event):
        # Update current_n label
        n = int(round(slider_n.get()))
        current_n.set(f"n={n}")
        
        # Update plots
        update_legendre()
    
    # Bind the sliders to the callback function
    slider_n.bind("<ButtonRelease-1>", slidermove_legendrepoly)
    
    # Initial plot
    update_legendre()
    
    # Start the Tkinter event loop
    root.mainloop()

def open_window_associated_legendre_poly():  

    
    # Create the root window for Tkinter
    root = tk.Tk()
    root.title("Associated Legendre Polynomial")
    root.state("zoomed")  # Start the window maximized
    
    # Parameters for the associated Legendre polynomial
    initial_m = 0  # Initial order
    initial_l = 0  # Initial degree
    
    # Create the plot figure
    fig = plt.figure(figsize=(8, 6))
    ax = fig.add_subplot(111, projection='3d')
    
    # Function to update the plot based on slider values
    def update_associated_poly(*args):
        m = int(slider_m.get())
        l = int(slider_l.get())
        
        # Update the labels next to the sliders
        label_value_m.config(text=f"{m}")
        label_value_l.config(text=f"{l}")
        
        if np.abs(m) > l:
            messagebox.showwarning("Warning", "Absolute Value of Order m cannot be greater than degree l!")
            root.update()  # Force the update of the messagebox window
            return  # Exit the function if the condition is invalid
    
        try:
            # Clear the current plot
            ax.clear()
    
            # Create a meshgrid for theta and phi
            theta = np.linspace(0, np.pi, 500)  # Polar angle (0 to pi)
            phi = np.linspace(0, 2 * np.pi, 500)  # Azimuthal angle (0 to 2pi)
            theta, phi = np.meshgrid(theta, phi)
    
            # Calculate the associated Legendre polynomial P^m_n(cos(theta))
            P_m_n = lpmv(m, l, np.cos(theta))
    
            # Convert spherical coordinates to Cartesian coordinates
            r = np.abs(P_m_n)  # Using the absolute value of the polynomial as the radial distance
            x = r * np.sin(theta) * np.cos(phi)
            y = r * np.sin(theta) * np.sin(phi)
            z = r * np.cos(theta)
    
            # Plotting the 3D surface
            ax.plot_surface(x, y, z, cmap='viridis', edgecolor='none')
    
            # Set plot labels and title
            ax.set_xlabel('X')
            ax.set_ylabel('Y')
            ax.set_zlabel('Z')
            ax.set_title(f'3D Plot of Associated Legendre Polynomial of m={m} and l={l}')
    
            # Redraw the canvas
            canvas.draw()
    
        except Exception as e:
            messagebox.showerror("Error", f"An error occurred while plotting: {e}")
            print(f"Error: {e}")
    
    # Frame to hold sliders and labels
    frame_sliders = ttk.Frame(root)
    frame_sliders.pack(side=tk.BOTTOM, fill=tk.X, padx=20, pady=10)
    
    # Label for m slider (Order)
    label_l = ttk.Label(frame_sliders, text="Degree (l):")
    label_l.pack(side=tk.LEFT, padx=(0, 10))
    
    # Slider for m (integer steps)
    slider_l = ttk.Scale(frame_sliders, from_=0, to=20, orient=tk.HORIZONTAL,length=500, command=update_associated_poly)
    slider_l.set(initial_m)
    slider_l.pack(side=tk.LEFT, padx=(0, 10))
    
    # Dynamic label to display selected value of m
    label_value_l = ttk.Label(frame_sliders, text=f"{initial_l}")
    label_value_l.pack(side=tk.LEFT)
    
    # Label for l slider (Degree)
    label_m = ttk.Label(frame_sliders, text="Order (m)")
    label_m.pack(side=tk.LEFT, padx=(20, 10))
    
    # Slider for l (integer steps)
    slider_m = ttk.Scale(frame_sliders, from_=-20, to=20, orient=tk.HORIZONTAL,length=500, command=update_associated_poly)
    slider_m.set(initial_l)
    slider_m.pack(side=tk.LEFT, padx=(0, 10))
    
    # Dynamic label to display selected value of l
    label_value_m = ttk.Label(frame_sliders, text=f"{initial_m}")
    label_value_m.pack(side=tk.LEFT)
    
    # Create a canvas to embed the plot in the Tkinter window
    canvas = FigureCanvasTkAgg(fig, master=root)
    canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=True)
    
    # Initial plot
    update_associated_poly()
    
    # Start the Tkinter event loop
    root.mainloop()
    #refer : https://www.mathworks.com/help/matlab/ref/legendre.html
    

def open_window_spherical_harmonics():
	# Create the root window for Tkinter
	root = tk.Tk()
	root.title("Spherical Harmonics Plot")
	root.state("zoomed")  # Start the window maximized
	
	# Parameters for the spherical harmonics
	initial_m = 0  # Initial order
	initial_l = 0  # Initial degree
	
	# Create the plot figure and the canvas for embedding it into Tkinter
	fig = plt.figure(figsize=(8, 6))
	ax = fig.add_subplot(111, projection='3d')
	canvas = FigureCanvasTkAgg(fig, master=root)
	
	# Function to update the plot based on slider values
	def update_spherical_harmonics(*args):
	    m = int(slider_m.get())  # Ensure that slider_m is defined before this call
	    l = int(slider_l.get())  # Ensure that slider_l is defined before this call
	
	    # Update the labels next to the s	liders	
	    label_value_m.config(text=f"{m}")	
	    label_value_l.config(text=f"{l}")	
	    	
	    if np.abs(m) > l:	
	        messagebox.showwarning("Warning", "Absolute value of Order m cannot be greater than 	degree l!")
	        root.update()  # Force the update of the messagebox window
	        return  # Exit the function if the condition is invalid
	
	    try:
	        # Clear the current plot
	        ax.clear()
	
	        # Create a meshgrid for theta and phi
	        theta = np.linspace(0, np.pi, 500)  # Polar angle (0 to pi)
	        phi = np.linspace(0, 2 * np.pi, 500)  # Azimuthal angle (0 to 2pi)
	        theta, phi = np.meshgrid(theta, phi)
	
	        # Calculate the spherical harmonics (Y_l^m(theta, phi))
	        Y_lm = sph_harm(m, l, phi, theta)
	
	        # Convert spherical coordinates to Cartesian coordinates
	        r = np.abs(Y_lm)  # Magnitude of the spherical harmonics
	        x = r * np.sin(theta) * np.cos(phi)
	        y = r * np.sin(theta) * np.sin(phi)
	        z = r * np.cos(theta)
	
	        # Plotting the 3D surface
	        ax.plot_surface(x, y, z, facecolors=plt.cm.viridis(r / r.max()), edgecolor='none')
	
	        # Set plot labels and title
	        ax.set_xlabel('X')
	        ax.set_ylabel('Y')
	        ax.set_zlabel('Z')
	        ax.set_title(f'Spherical Harmonics for m={m} and l={l}')
	
	        # Redraw the canvas
	        canvas.draw()
	
	    except Exception as e:
	        messagebox.showerror("Error", f"An error occurred while plotting: {e}")
	        print(f"Error: {e}")
	
	# Frame to hold sliders and labels
	frame_sliders = ttk.Frame(root)
	frame_sliders.pack(side=tk.BOTTOM, fill=tk.X, padx=20, pady=10)
	
	# Label for l slider (Degree)
	label_l = ttk.Label(frame_sliders, text="Degree (l):")
	label_l.pack(side=tk.LEFT, padx=(0, 10))
	
	# Slider for l (integer steps)	
	slider_l = ttk.Scale(frame_sliders, from_=0, to=10, orient=tk.HORIZONTAL, length=400, 	command=update_spherical_harmonics)
	slider_l.set(initial_l)
	slider_l.pack(side=tk.LEFT, padx=(0, 10))
	
	# Dynamic label to display selected value of l
	label_value_l = ttk.Label(frame_sliders, text=f"{initial_l}")
	label_value_l.pack(side=tk.LEFT)
	
	# Label for m slider (Order)
	label_m = ttk.Label(frame_sliders, text="Order (m):")
	label_m.pack(side=tk.LEFT, padx=(20, 10))
	
	# Slider for m (integer steps)	
	slider_m = ttk.Scale(frame_sliders, from_=-10, to=10, orient=tk.HORIZONTAL, length=400, 	command=update_spherical_harmonics)
	slider_m.set(initial_m)
	slider_m.pack(side=tk.LEFT, padx=(0, 10))
	
	# Dynamic label to display selected value of m
	label_value_m = ttk.Label(frame_sliders, text=f"{initial_m}")
	label_value_m.pack(side=tk.LEFT)
	
	# Embed the plot canvas into the Tkinter window
	canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=True)
	
	# Initial plot
	update_spherical_harmonics()
	
	# Start the Tkinter event loop
	root.mainloop()

def open_window_tunneling():
    hbar = 1  # reduced Planck's constant, J·s
    me = 1    # mass of an electron, kg
    
    N = 1000  # Number of points for plotting
    initial_a = 1.0
    initial_E = 0
    initial_V = 10.0
    message_shown = False  # To track if the message has been shown

    def calculate_wavefunctions(a, E, V):
        k = np.sqrt(2 * me * E) / hbar

        # Ensure alpha is real for E < V, otherwise handle complex case
        if E <= V:
            alpha = np.sqrt(2 * me * (V - E)) / hbar
        else:
            alpha = np.sqrt(2 * me * np.abs(V - E)) * 1j / hbar

        # Coefficients for wavefunction inside and outside the barrier
        A=1
        F = (2 * k * alpha * np.exp(-1j * k * a))*A / ((alpha**2 - k**2) * np.cosh(alpha * a) + 2 * 1j * k * alpha * np.sinh(alpha * a))
        D = (F / (2 * alpha)) * (alpha + 1j * k) * np.exp(1j * k * a - alpha * a)
        C = (F / (2 * alpha)) * (alpha - 1j * k) * np.exp(1j * k * a + alpha * a)
        B = C + D - A

        # Combine regions into a single x array
        x1 = np.linspace(-10, 0, N)
        x2 = np.linspace(0, a, N)
        x3 = np.linspace(a, 10, N)

        x = np.concatenate((x1, x2, x3))  # Combined x array

        if V == 0 or a == 0:
            # When V=0 or a=0, it's a single continuous wave
            psi = np.exp(1j * k * x)  # Same wave across the entire region
        else:
            # Calculate wavefunctions for different regions
            psi_1 = A*np.exp(1j * k * x1) + B * np.exp(-1j * k * x1)
            psi_2 = C * np.exp(-alpha * x2) + D * np.exp(alpha * x2)
            psi_3 = F * np.exp(1j * k * x3)
            
            # Combine wavefunctions into a single array
            psi = np.concatenate((psi_1, psi_2, psi_3))

        return x, psi

    def update_Tunneling():
        nonlocal message_shown  # Allow the function to modify the flag

        a = float(slider_a.get())
        E = float(slider_E.get())
        V = float(slider_V.get())

        # Check if either 'a' or 'V' is zero and automatically set the other to zero
        if a == 0 or V == 0:
            if a != 0:
                slider_a.set(0)
                a = 0
            if V != 0:
                slider_V.set(0)
                V = 0
            if not message_shown:  # Only show the message once
                messagebox.showinfo("Barrier Vanishes", "When either width (a) or potential (V) is zero, the barrier vanishes.")
                message_shown = True

        # Get the wavefunctions
        x, psi = calculate_wavefunctions(a, E, V)

        k = np.sqrt(2 * me * E / hbar**2)
        
        # Ensure alpha is real for E < V, otherwise handle complex case
        alpha = np.sqrt(2 * me * np.abs(V - E) / hbar**2) * (1 if E < V else 1j)

        # Special case for V = 0
        if V == 0:
            T = 1.0
            R = 0.0
        else:
            T = (4 * k**2 * alpha**2) / ((alpha**2 - k**2) * (np.sinh(alpha * a))**2 + 4 * k**2 * alpha**2 * (np.cosh(alpha * a))**2)
            R = 1 - T
        
        # Clear previous plots
        plot_tunnel.clear()
        
        current_a.set(f"a = {a:.2f}")
        current_E.set(f"E = {E:.2f}")
        current_V.set(f"V = {V:.2f}")
        plot_tunnel.axhline(y=0, color='black', linestyle='--')  # y=0
        
        
        

        # Plot the wave function (imaginary part)
        plot_tunnel.plot(x, np.imag(psi), label=r"Im($\psi(x)$)", color="b")

        # Plot the reflection and transmission as constant lines
        plot_tunnel.axhline(y=np.real(R), linestyle='--', color='cyan', label='R (Reflection Coefficient)')
        plot_tunnel.axhline(y=np.real(T), linestyle='--', color='darkblue', label='T (Transmission Coefficient)')
        
        
        # Add text annotations for R and T
        plot_tunnel.text(0, np.real(R), f'R = {np.real(R):.4f}', fontsize=10, color='cyan', ha='left')
        plot_tunnel.text(0, np.real(T), f'T = {np.real(T):.4f}', fontsize=10, color='darkblue', ha='left')
        
        # Plot the width and potential barrier lines
        if a > 0 and V > 0:
            # Horizontal line for potential
            plot_tunnel.plot([0, a], [V, V], color='black', linewidth=2, label='Potential Barrier')
            # Vertical lines for width
            plot_tunnel.plot([0, 0], [0, V], color='black', linewidth=2)  # Start of the barrier
            plot_tunnel.plot([a, a], [0, V], color='black', linewidth=2)  # End of the barrier
            
            plot_tunnel.plot([0, a], [-V, -V], color='g', linewidth=2, linestyle='--', label='Magnitude of Potential Barrier')
            # Vertical lines for width
            plot_tunnel.plot([0, 0], [0, -V], color='g', linewidth=2,linestyle='--')  # Start of the barrier
            plot_tunnel.plot([a, a], [0, -V], color='g', linewidth=2,linestyle='--')  # End of the barrier

        # Set titles and labels
        plot_tunnel.set_title(f'(Tunnel Width = {a:.2f}, Energy = {E:.2f}, Potential = {V:.2f})')
        plot_tunnel.set_xlabel('Position x (arbitrary units)')
        plot_tunnel.set_ylabel(r"$\psi(x)$")
        plot_tunnel.legend()
        plot_tunnel.grid(True)
        
        # Fix y-axis limits based on max(V)
        max_y = max(10, V + 2)  # Ensure there's a minimum range, so it doesn't squeeze at low V
        plot_tunnel.set_ylim(-max_y, max_y)  # Symmetric limits around 0

        # Set y-axis ticks with steps of 1
        y_ticks = np.arange(-max_y, max_y+1, 1)  # y-ticks with step of 1
        plot_tunnel.set_yticks(y_ticks)

        canvas.draw()

    # Automatically set the other slider to 1 if it's 0 and the user moves one
    def on_slider_change_a(value):
        a = float(slider_a.get())
        V = float(slider_V.get())
        if V == 0 and a > 0:
            slider_V.set(1)  # Automatically set V to 1
        update_Tunneling()

    def on_slider_change_V(value):
        a = float(slider_a.get())
        V = float(slider_V.get())
        if a == 0 and V > 0:
            slider_a.set(1)  # Automatically set a to 1
        update_Tunneling()
    
    # Create the main window
    root = tk.Tk()
    root.title("Quantum Tunneling")
    root.state('zoomed')
    
    # Create a single figure for both plots
    fig, plot_tunnel = plt.subplots(figsize=(8, 6))
    
    # Embed the figure into the Tkinter window
    canvas = FigureCanvasTkAgg(fig, master=root)
    canvas.draw()
    canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=1)
    
    # Frame to hold sliders and labels
    frame_sliders = ttk.Frame(root)
    frame_sliders.pack(side=tk.TOP, fill=tk.X, padx=20, pady=10)
    
    # Slider for constant 'a'
    label_a = ttk.Label(frame_sliders, text="Width (a):")
    label_a.pack(side=tk.LEFT, padx=(0, 10))
    slider_a = ttk.Scale(frame_sliders, from_=0, to=10.0, orient=tk.HORIZONTAL, length=300, command=on_slider_change_a)
    slider_a.set(initial_a)
    slider_a.pack(side=tk.LEFT, padx=(0, 10))
    
    current_a = tk.StringVar()
    current_a.set(f"a = {initial_a:.2f}")
    label_current_a = ttk.Label(frame_sliders, textvariable=current_a)
    label_current_a.pack(side=tk.LEFT, padx=(0, 20))
    
    # Slider for energy 'E'
    label_b = ttk.Label(frame_sliders, text="Energy (E):")
    label_b.pack(side=tk.LEFT, padx=(0, 10))
    slider_E = ttk.Scale(frame_sliders, from_=0, to=10.0, orient=tk.HORIZONTAL, length=300, command=lambda val: update_Tunneling())
    slider_E.set(initial_E)
    slider_E.pack(side=tk.LEFT, padx=(0, 10))
    
    current_E = tk.StringVar()
    current_E.set(f"E = {initial_E:.2f}")
    label_current_E = ttk.Label(frame_sliders, textvariable=current_E)
    label_current_E.pack(side=tk.LEFT, padx=(0, 20))
    
    # Slider for potential 'V'
    label_c = ttk.Label(frame_sliders, text="Potential (V):")
    label_c.pack(side=tk.LEFT, padx=(0, 10))
    slider_V = ttk.Scale(frame_sliders, from_=0, to=20.0, orient=tk.HORIZONTAL, length=300, command=on_slider_change_V)
    slider_V.set(initial_V)
    slider_V.pack(side=tk.LEFT, padx=(0, 10))
    
    current_V = tk.StringVar()
    current_V.set(f"V = {initial_V:.2f}")
    label_current_V = ttk.Label(frame_sliders, textvariable=current_V)
    label_current_V.pack(side=tk.LEFT, padx=(0, 20))
    
    update_Tunneling()  # Initial plot update
    root.mainloop()

import tkinter as tk
from tkinter import ttk
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg

def open_window_free_electron_gas():
    # Constants
    hbar = 1  # Reduced Planck's constant
    me = 1    # Mass of an electron
    initial_q = 1
    initial_N = 1
    initial_ymax = 20
    y_ultimate = 10000
    initial_x = 10
    x_ultimate = 1000

    # Calculate energy and pressure as functions of V, with q and N as parameters
    def calculate_Energy(V, q, N):
        constant_fac = (hbar**2) * ((3 * np.pi**2)**(5 / 3)) / (10 * me * np.pi**2)
        variable_fac = (N * q)**(5 / 3)
        E = constant_fac * variable_fac * V**(-2 / 3)
        return E

    def calculate_Pressure(V, q, N):
        rho = N * q  # Density
        P = (3 * np.pi**2)**(2 / 3) * (hbar**2) * (1 / (5 * me)) * rho**(5 / 3) * V**(-5 / 3)
        return P

    def update_free_electron_gas():
        # Round values to ensure integer steps
        q = round(float(slider_q.get()))
        N = round(float(slider_N.get()))
        ymax = round(float(slider_ymax.get()))
        xmax = round(float(slider_xmax.get()))

        # Generate a range for volume V to plot E and P against, using xmax from slider
        V_values = np.linspace(0.001, xmax, 10000)
        E_values = calculate_Energy(V_values, q, N)
        P_values = calculate_Pressure(V_values, q, N)

        # Clear previous plots
        plot_free_electron.clear()

        current_q.set(f"q = {q}")
        current_N.set(f"N = {N}")
        current_y.set(f"{ymax}")  # Update current_y with ymax
        current_x.set(f"{xmax}")  # Update current_x with xmax

        # Plot energy and pressure as functions of V
        plot_free_electron.plot(V_values, E_values, color='cyan', label='E (Energy)')
        plot_free_electron.plot(V_values, P_values, color='red', label='P (Pressure)')

        # Set titles and labels
        plot_free_electron.set_title(f'Free Electron Gas (q={q}, N={N})')
        plot_free_electron.set_xlabel('Volume (V)')
        plot_free_electron.set_ylabel('Energy and Pressure')
        plot_free_electron.legend()
        plot_free_electron.grid(True)

        # Adjust y-limits for consistent scale
        plot_free_electron.set_ylim(0, ymax)
        
        canvas.draw()

    # Tkinter window setup
    root = tk.Tk()
    root.title("Free Electron Gas")
    root.state('zoomed')

    # Matplotlib figure setup
    fig, plot_free_electron = plt.subplots(figsize=(8, 6))
    plt.tight_layout(pad=1.0)
    canvas = FigureCanvasTkAgg(fig, master=root)
    canvas.draw()
    canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=1)

    # Frame to hold sliders in a grid format
    frame_sliders = ttk.Frame(root)
    frame_sliders.pack(side=tk.TOP, fill=tk.X, padx=5, pady=5)

    # First row, first column: Charge Density (q)
    frame_q = ttk.Frame(frame_sliders)
    frame_q.grid(row=0, column=0, padx=5, pady=5)
    ttk.Label(frame_q, text="Charge Density (q):").grid(row=0, column=0, sticky="w")
    slider_q = ttk.Scale(frame_q, from_=1, to=10, orient=tk.HORIZONTAL, length=500, command=lambda val: update_free_electron_gas())
    slider_q.set(initial_q)
    slider_q.grid(row=0, column=1, sticky="w")
    current_q = tk.StringVar()
    current_q.set(f"q = {initial_q}")
    ttk.Label(frame_q, textvariable=current_q).grid(row=0, column=2, sticky="w", padx=10)

    # First row, second column: Number of Particles (N)
    frame_N = ttk.Frame(frame_sliders)
    frame_N.grid(row=0, column=1, padx=5, pady=5)
    ttk.Label(frame_N, text="Number of Particles (N):").grid(row=0, column=0, sticky="w")
    slider_N = ttk.Scale(frame_N, from_=1, to=100, orient=tk.HORIZONTAL, length=500, command=lambda val: update_free_electron_gas())
    slider_N.set(initial_N)
    slider_N.grid(row=0, column=1, sticky="w")
    current_N = tk.StringVar()
    current_N.set(f"N = {initial_N}")
    ttk.Label(frame_N, textvariable=current_N).grid(row=0, column=2, sticky="w", padx=10)

    # Second row, first column: Maximum y
    frame_y = ttk.Frame(frame_sliders)
    frame_y.grid(row=1, column=0, padx=5, pady=5)
    ttk.Label(frame_y, text="Maximum y:").grid(row=0, column=0, sticky="w")
    slider_ymax = ttk.Scale(frame_y, from_=initial_ymax, to=y_ultimate, orient=tk.HORIZONTAL, length=500, command=lambda val: update_free_electron_gas())
    slider_ymax.set(initial_ymax)
    slider_ymax.grid(row=0, column=1, sticky="w")
    current_y = tk.StringVar()
    current_y.set(f"{initial_ymax}")
    ttk.Label(frame_y, textvariable=current_y).grid(row=0, column=2, sticky="w", padx=10)

    # Second row, second column: Maximum x
    frame_x = ttk.Frame(frame_sliders)
    frame_x.grid(row=1, column=1, padx=5, pady=5)
    ttk.Label(frame_x, text="Maximum x:").grid(row=0, column=0, sticky="w")
    slider_xmax = ttk.Scale(frame_x, from_=initial_x, to=x_ultimate, orient=tk.HORIZONTAL, length=500, command=lambda val: update_free_electron_gas())
    slider_xmax.set(initial_x)
    slider_xmax.grid(row=0, column=1, sticky="w")
    current_x = tk.StringVar()
    current_x.set(f"{initial_x}")
    ttk.Label(frame_x, textvariable=current_x).grid(row=0, column=2, sticky="w", padx=10)

    # Initial plot update
    update_free_electron_gas()
    root.mainloop()
    
#-------------
def open_window_distribution_functions():
    # Constants
    k_B = 1  # Boltzmann constant in arbitrary units
    initial_alpha = 1  # Initial value for alpha (scaling factor for mu)
    initial_T = 3  # Initial temperature in Kelvin

    # Functions to calculate distributions
    def maxwell_boltzmann(E, mu, T):
        return np.exp(-(E - mu) / (k_B * T))

    def fermi_dirac(E, mu, T):
        return 1 / (np.exp((E - mu) / (k_B * T)) + 1)

    def bose_einstein(E, mu, T):
        return np.where(E > mu, 1 / (np.exp((E - mu) / (k_B * T)) - 1), np.inf)

    def update_plot(*args):
        # Get values from sliders
        T = slider_T.get()
        alpha = slider_alpha.get()
        mu = alpha * k_B * T
        E = np.linspace(-5, 50, 500)  # Energy range

        # Update labels
        current_T.set(f"T = {T:.2f} K")
        current_alpha.set(f"α = {alpha:.2f}")

        # Clear the previous plot
        ax.clear()

        # Plot distributions
        ax.plot(E, maxwell_boltzmann(E, mu, T), label="Maxwell-Boltzmann", color="blue")
        ax.plot(E, fermi_dirac(E, mu, T), label="Fermi-Dirac", color="red")
        ax.plot(E, bose_einstein(E, mu, T), label="Bose-Einstein", color="green")

        # Labels and title
        ax.set_xlabel("Energy (E)")
        ax.set_ylabel("Distribution Function n(E)")
        ax.set_title(f"Distributions at T = {T:.2f} K and μ = {mu:.2f}")
        ax.legend()
        ax.grid(True)

        # Set y-axis limits
        ax.set_ylim(0, 2)

        # Redraw the canvas
        canvas.draw()

    # Create the main window
    root = tk.Tk()
    root.title("Distribution Functions: Maxwell-Boltzmann, Fermi-Dirac, Bose-Einstein")
    root.state('zoomed')

    # Create a figure and axis for the plot
    fig, ax = plt.subplots(figsize=(8, 6))
    plt.tight_layout(pad=1.0)

    # Embed the plot into the Tkinter window
    canvas = FigureCanvasTkAgg(fig, master=root)
    canvas.draw()
    canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=1)

    # Frame to hold sliders in a grid format
    frame_sliders = ttk.Frame(root)
    frame_sliders.pack(side=tk.TOP, fill=tk.X, padx=5, pady=5)

    # Temperature (T) slider
    frame_T = ttk.Frame(frame_sliders)
    frame_T.grid(row=0, column=0, padx=5, pady=5)
    ttk.Label(frame_T, text="Temperature (T) [K]:").grid(row=0, column=0, sticky="w")
    slider_T = ttk.Scale(frame_T, from_=0, to=10.0, orient=tk.HORIZONTAL, length=550, command=update_plot)
    slider_T.set(initial_T)
    slider_T.grid(row=0, column=1, sticky="w")
    current_T = tk.StringVar()
    current_T.set(f"T = {initial_T:.2f} K")
    ttk.Label(frame_T, textvariable=current_T).grid(row=0, column=2, sticky="w", padx=10)

    # Alpha (α) slider
    frame_alpha = ttk.Frame(frame_sliders)
    frame_alpha.grid(row=0, column=1, padx=5, pady=5)
    ttk.Label(frame_alpha, text="Alpha (α):").grid(row=0, column=0, sticky="w")
    slider_alpha = ttk.Scale(frame_alpha, from_=0, to=3.0, orient=tk.HORIZONTAL, length=550, command=update_plot)
    slider_alpha.set(initial_alpha)
    slider_alpha.grid(row=0, column=1, sticky="w")
    current_alpha = tk.StringVar()
    current_alpha.set(f"α = {initial_alpha:.2f}")
    ttk.Label(frame_alpha, textvariable=current_alpha).grid(row=0, column=2, sticky="w", padx=10)

    # Initial plot update
    update_plot()

    # Start the main loop
    root.mainloop()




##########################################################################################
    
# Main application window
root = tk.Tk()
root.title("Quantum Animations")
root.state('zoomed')

# Create a frame to hold the buttons
button_frame = ttk.Frame(root)
button_frame.pack(pady=20)

# List of button configurations: (row, column, text, command)
buttons = [
    (0, 0, "Infinite Square Well", open_window_infinitesquarewell),
    (0, 1, "Harmonic Oscillator", open_window_HO),
    (0, 2, "Delta Dirac Potential Well Bound State", open_window_delta_dirac_2D),
    (0, 3, "Delta Dirac Potential Well Scattering State", open_window_delta_scatter),
    (0, 4, "Legendre Polynomials", open_window_legendre_poly),
    (1, 0, "Associated Legendre Polynomials", open_window_associated_legendre_poly),
    (1, 1, "Spherical Harmonics", open_window_spherical_harmonics),
    (1, 2, "Quantum Tunneling", open_window_tunneling),
    (1, 3, "Free Electron Gas", open_window_free_electron_gas),
    (1, 4, "Distributions", open_window_distribution_functions)
]

# Create and place buttons in the grid
for row, col, text, command in buttons:
    button = ttk.Button(button_frame, text=text, width=40, command=command)
    button.grid(row=row, column=col, padx=10, pady=10, ipadx=20, ipady=10, sticky="nsew")

# Configure grid to ensure uniform button sizes
for i in range(5):  # 5 columns
    button_frame.grid_columnconfigure(i, weight=1)
for i in range(2):  # 2 rows
    button_frame.grid_rowconfigure(i, weight=1)

# Start the tkinter event loop
root.mainloop()

