# DensityCube

`DensityCube` is a Python package for working with 3D density cube files. This package allows you to load, process, and analyze density grids, as well as calculate various properties like gradients and integrated densities. The code is tested using unit tests to ensure the correctness of the functionality.

## Features
- Load and parse cube files
- Access voxel dimensions and grid shape
- Calculate derivatives and Hessians for 3D grids
- Integrate density over the grid
- Unit tests to verify functionality

### MonteCarlo
`MonteCarlo` provides functionality to perform Monte Carlo integration simulations and visualize the results. It is a versatile tool for numerical computations requiring random sampling.

#### Key features:
- Perform Monte Carlo integration over a user-specified range and step size.
- Save integration results to temporary files for quick analysis.
- Generate and save plots of integration results to visualize convergence.
- Fully tested with integration results verified through unit tests.

#### Usage:
1. **Command-line Execution:**
   Execute Monte Carlo integration and generate a plot:
   ```bash
   python monteCarlo.py --xrange 0 1 --yrange 0 1 --nsteps 1000 --stepsize 100
   ```


## Installation
1. Clone the repository:
```bash
   git clone https://github.com/jarcoaleksandra/numerical.git
   cd numerical
```

2. Create a virtual environment:
```bash
   python3 -m venv venv
```

3. Activate the virtual environment:

   -For Unix/macOS:
   ```bash
   source venv/bin/activate
   ```	
   	
   -For Windows:
   ```bash
   venv\Scripts\activate
   ```
4. Install required dependencies:
	```bash
	pip install -r requirements.txt
	```
## Usage
### Example for DensityCube
```bash
cd example
python example_script.py
```

### Example for MonteCarlo
```bash
python monteCarlo.py --xrange 0 1 --yrange 0 1 --nsteps 1000 --stepsize 100
```

