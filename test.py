import tkinter as tk
from tkinter import ttk
import numpy as np
import matplotlib.pyplot as plt
import tkinter as tk
import math
from scipy.special import hermite
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg

# Function to open a new window for Button A
def open_window_a():
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
    
 

# Function to open a new window for Button B
def open_window_b():
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
    
def open_window_delta_dirac_2D():
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
    slider_n = ttk.Scale(frame_sliders, from_=0, to=100, orient=tk.HORIZONTAL, command=lambda val:     update_Dirac())    
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
    
    
    
    
    
# Main application window
root = tk.Tk()
root.title("Main Window")
root.geometry("300x200")

# Create buttons for A and B
button_a = ttk.Button(root, text="Button A", command=open_window_a)
button_a.pack(pady=20)

button_b = ttk.Button(root, text="Button B", command=open_window_b)
button_b.pack(pady=20)

button_c = ttk.Button(root, text="Delta Dirac Potential Well 2D", command=open_window_delta_dirac_2D)
button_c.pack(pady=20)
# Start the tkinter event loop
root.mainloop()
