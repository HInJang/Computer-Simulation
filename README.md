# Computer-Simulation

A Python-based orbital simulation of the solar system.
This project uses numerical integration methods (Beeman, Euler, and Symplectic Euler) to model the orbits, total energy, and planetary alignments of celestial bodies in a simplified solar system.

### 📁 Project Structure
- `planet.py`: Defines the Planet class and handles position, velocity, and acceleration updates.

- `solar.py`: Runs the main simulation logic, visualization, energy calculations, and analysis.

- `parameters-solar.txt`: Input data file containing simulation parameters and solar system data.

### 🚀 Features
- **Orbital Simulation**: Animates planetary orbits using matplotlib.

- **Energy Analysis**: Compares total energy conservation between Beeman and Euler-based methods.

- **Planetary Alignment Detection**: Detects alignment based on angular thresholds.

- **Orbital Periods**: Compares simulated vs. actual orbital periods using regression analysis.

### 🧪 How to Run

<pre> python solar.py </pre>
- Upon execution, the simulation will:

- Display animated planetary orbits

- Plot orbital periods and regression fit

- Visualize energy conservation over time (Beeman vs. Euler)

- Detect and display planetary alignments

### 📘 Numerical Methods Used
- Beeman Integration

- Direct Euler Method

- Symplectic Euler Method

### ⚙️ Requirements
Install dependencies with:

<pre> pip install numpy matplotlib seaborn scikit-learn  </pre>

### 📊 Data Source
Planetary data taken from: NASA JPL Horizons System
*(see `parameters-solar.txt` for details)*

### ⚠️ Notes
Units: AU (astronomical units), M⊕ (Earth mass), yr (Earth year)

Conversion factors included for energy output in SI units (Joules)
