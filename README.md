# Origins of Instability in Dynamical Systems on Undirected Networks

## Project Overview

This project analyzes the origins of instability in dynamical systems on undirected networks. It investigates the relationship between network topology, dynamical parameters, and system stability through both analytical and numerical approaches. The work focuses on SIS (Susceptible-Infected-Susceptible) epidemic models and their stability properties on various network structures.

---

## Core Symbols and Parameters

### Network Parameters

| Symbol | Description | Type | Example Value |
|--------|-------------|------|---------------|
| `N` | Number of nodes in the network | int | 200, 500, 5000 |
| `n` | Number of nodes (alternative notation) | int | 500 |
| `p` | Probability for Erdős-Rényi graph generation | float | 0.02 |
| `m` | Number of edges to attach per node in Barabási-Albert model | int | 5 |
| `G` | NetworkX Graph object | nx.Graph | - |
| `A` | Adjacency matrix | np.ndarray | 2D array |
| `degree` | Node degree vector | np.ndarray | 1D array |
| `avg` | Average degree of network | float | Computed from degrees |
| `num_nodes` | Number of nodes in subgraph | int | - |
| `num_edges` | Number of edges in graph | int | - |
| `s` | Sparsity/density parameter | float | num_edges / total_possible_edges |

### Dynamical System Parameters

| Symbol | Description | Type | Example Value |
|--------|-------------|------|---------------|
| `beta` (β) | Infection transmission rate | float | 0.0 - 0.5 |
| `alpha` (α) | Recovery rate | float | 1 |
| `chi` (χ) | Coupling parameter | float | -0.1, 0.01, 0.3 |
| `mu` (μ) | Power-law exponent for degree dependence | float | 1 |
| `nu` (ν) | Power-law exponent (off-diagonal) | float | -1 |
| `rho` (ρ) | Power-law exponent (node interaction) | float | 0 |
| `mu_net` | Network coupling constant | float | 1 |

### State Variables

| Symbol | Description | Type | Units |
|--------|-------------|------|-------|
| `x` | Infection state variable | np.ndarray | [0, 1] |
| `x_avg` | Average/mean infection state | float | Steady-state value |
| `x_star` | Fixed point / equilibrium value | float | Steady-state solution |
| `x_final` | Final state after integration | np.ndarray | - |
| `x_0` | Initial condition | np.ndarray | Random or specified |
| `x0` | Initial state vector | np.ndarray | Random [0, 1] |

### Eigenvalue and Stability Symbols

| Symbol | Description | Type | Notes |
|--------|-------------|------|-------|
| `lambda` (λ) | Eigenvalue | complex | Real part determines stability |
| `lambda_max` | Maximum eigenvalue | float | Largest eigenvalue |
| `lam_1` | Largest eigenvalue of M | list | Analytical |
| `lam_11` | Largest eigenvalue of M1 | list | Linearized system |
| `lam_t` | Theoretical eigenvalue (mean approximation) | list | Derived from theory |
| `lam_t1` | Theoretical eigenvalue (linearized) | list | Enhanced prediction |
| `lam_P1` | Eigenvalue prediction (max degree) | list | λ_max ≈ -α - βx*k_max |
| `lam_P2` | Eigenvalue prediction (min degree) | list | λ_min ≈ -α - βx*k_min |
| `lam_numerical` | Numerically computed largest eigenvalue | float | From Jacobian at steady state |

### Jacobian and Stability Matrices

| Symbol | Description | Dimensions | Purpose |
|--------|-------------|-----------|---------|
| `M` | Linearized dynamics matrix (mean-field) | (n, n) | M_ij = β·A_ij, M_ii = -α |
| `M1` | Linearized dynamics matrix (linearized equilibrium) | (n, n) | Corrected for heterogeneous state |
| `J` | Jacobian matrix at steady state | (n, n) | J_ij = ∂(dx_i/dt)/∂x_j |
| `values` | Eigenvalues of M | 1D array | Real parts from eig(M) |
| `values1` | Eigenvalues of M1 | 1D array | Real parts from eig(M1) |
| `eigvals` | Eigenvalues from eigsh | 1D array | Sparse matrix eigenvalues |
| `vect` | Eigenvectors from dense eig | (n, n) | Columns are eigenvectors |
| `eigvecs` | Eigenvectors from sparse eigsh | (n, k) | Sparse eigenvectors |

---

## Key Mathematical Models

### SIS Dynamical Equation

```
dx_i/dt = -x_i + β(1 - x_i)(A @ x)_i
```

**Key Components:**
- `dxdt(t, x, beta, A)` - Function definition for SIS dynamics
- `x[i]` - Infection probability of node i
- `A @ x` - Convolution of adjacency matrix with infection state

### Fixed Point Analysis

At equilibrium (dx/dt = 0):

```
x_avg = 1 - α/(β·k_avg)
x_star = 1 - 1/(1 + (β - 1/k_avg)·k_avg)
```

---

## Core Functions

### Utility Functions

#### `calculate_mean_non_zero_off_diagonal(matrix)` 
- **Location:** [dynimics_fx_pt_stability.py](dynimics_fx_pt_stability.py), [fig_3.py](fig_3.py)
- **Parameters:** `matrix` (np.ndarray) - Input matrix
- **Returns:** float - Mean of all non-zero off-diagonal elements
- **Purpose:** Compute mean coupling strength for stability analysis

#### `calculate_degree_sparsity(G)`
- **Location:** [dynimics_fx_pt_stability.py](dynimics_fx_pt_stability.py), [fig_3.py](fig_3.py)
- **Parameters:** `G` (nx.Graph) - NetworkX graph object
- **Returns:** float - Sparsity metric = ⟨k²⟩/⟨k⟩²
- **Purpose:** Quantify network heterogeneity

### Dynamics Functions

#### `dxdt(t, x, beta, A)`
- **Location:** [dynimics_fx_pt_stability.py](dynimics_fx_pt_stability.py), [SIS_FX_PT.py](SIS_FX_PT.py)
- **Parameters:** 
  - `t` (float) - Time
  - `x` (np.ndarray) - State vector
  - `beta` (float) - Infection rate
  - `A` (np.ndarray) - Adjacency matrix
- **Returns:** np.ndarray - Time derivative dx/dt
- **Equation:** `-x + β(1 - x) ⊙ (A @ x)`

---

## Data Files

### Network Data Files

| File | Format | Description | Dimensions |
|------|--------|-------------|-----------|
| `A_20_35.mat` | MATLAB | Adjacency matrix | 20×20 or 35×35 |
| `A_500_avg30.mat` | MATLAB | Adjacency matrix with avg degree ≈ 30 | 500×500 |
| `A_BA_avgdeg12.mat` | MATLAB | Barabási-Albert network, avg degree = 12 | 200×200 |
| `A_BA_avgdeg6.mat` | MATLAB | Barabási-Albert network, avg degree = 6 | 200×200 |
| `BA_m_5_n_200.txt` | Text | Barabási-Albert (m=5, n=200) adjacency matrix | 200×200 |

### Results Data Files

| File | Format | Columns | Description |
|------|--------|---------|-------------|
| `data_sis_beta_3_chi_1_change_tim_av_10_5000_0.5_2_txt` | Text | n, s, β, χ, λ_1, λ_out, λ_out1, λ_out2 | SIS simulation results |

---

## File Descriptions

### [dynimics_fx_pt_stability.py](dynimics_fx_pt_stability.py)
**Purpose:** Analyze fixed point stability through eigenvalue computation

**Key Variables:**
- `mu_net`, `chi`, `alpha` - System parameters
- `n`, `p` - Network generation parameters
- `data` - β values array (0.01 to 0.5)
- `lam_t`, `lam_1` - Theoretical vs. numerical eigenvalues
- `lam_P1`, `lam_P2` - Eigenvalue predictions from max/min degree

**Main Operations:**
1. Generate Erdős-Rényi network
2. Extract giant connected component
3. For each β value:
   - Construct linearization matrices M, M1
   - Compute eigenvalues analytically
   - Integrate dynamics to steady state
   - Compute Jacobian eigenvalues numerically
4. Plot eigenvalue comparison

---

### [SIS_FX_PT.py](SIS_FX_PT.py)
**Purpose:** Simulate SIS dynamics and validate theoretical predictions

**Key Variables:**
- `beta_values` - Range [0, 0.5]
- `T_max` - Maximum integration time (10000)
- `t_eval` - Time evaluation points (10000 steps)
- `N` - Number of nodes (500)
- `theo_avg`, `theo_star` - Theoretical equilibrium predictions
- `x_star_values` - Numerically computed equilibrium values

**Main Operations:**
1. Load network from MATLAB file
2. Compute average degree and theoretical equilibrium
3. For each β:
   - Random initial condition
   - Integrate dynamics using solve_ivp
   - Extract final state and compute mean
4. Compare theoretical vs. numerical equilibrium

---

### [fig_1.py](fig_1.py)
**Purpose:** Create phase space plot identifying bifurcation transitions

**Key Variables:**
- `red_points`, `blue_points`, `green_points` - Point sets by category
- `g_zero` - Critical χ value where eigenvalue crosses zero
- `lam_1_vals` - Eigenvalue data
- `lambda_out` - Mean-field prediction
- `lambda_out1`, `lambda_out2` - Alternative predictions
- `selected_points` - Dictionary of one point per (β, n) pair

**Classification Logic:**
- **Red:** λ_numerical ≈ λ_mean
- **Blue:** λ_numerical ≈ λ_max_degree
- **Green:** λ_numerical ≈ λ_min_degree

---

### [fig_2(1).py](fig_2(1).py)
**Purpose:** Scatter plot of eigenvector components vs. node degree

**Key Variables:**
- `J` - Jacobian matrix with heterogeneous structure
- `eigvals`, `eigvecs` - Eigenvalues and eigenvectors from eigsh
- `v` - Absolute values of largest eigenvector components
- `theo` - Theoretical prediction: k^(ν+1)
- `chi_values` - Three chi values for comparison

**Matrix Structure:**
```
J[i,j] = (k_i^ν)(k_j^ρ)A[i,j]
J[i,i] = -β - χ·k_i^μ
```

---

### [fig_2(2).py](fig_2(2).py)
**Purpose:** Network visualization colored by eigenvector values

**Key Variables:**
- `H` - Induced subgraph (nodes with degree > median)
- `high_degree_nodes` - Filtered node set
- `degree_threshold` - 50th percentile of degree distribution
- `node_colors` - Eigenvector component values
- `node_sizes` - Proportional to degree

**Visualization:**
- Color map: 'coolwarm' (blue = low, red = high eigenvector values)
- Node size proportional to degree
- Black edges showing connectivity

---

### [fig_3.py](fig_3.py)
**Purpose:** Extended analysis of fixed point stability with multiple parameters

**Key Variables:**
- Similar to [dynimics_fx_pt_stability.py](dynimics_fx_pt_stability.py)
- Generates comprehensive eigenvalue landscape

---

## Theoretical Relationships

### Stability Conditions

The system is stable when λ_max < 0:

```
λ_max < 0  ⟹  System approaches equilibrium
λ_max > 0  ⟹  System is unstable (bifurcation)
λ_max = 0  ⟹  Critical transition point
```

### Eigenvalue Approximations

1. **Mean-field approximation:**
   ```
   λ_mean = -α + β·⟨A_off-diag⟩·(n-1)·s·sparsity
   ```

2. **Degree-dependent approximation:**
   ```
   λ_max_deg = -α - β·x*·max(degree)
   λ_min_deg = -α - β·x*·min(degree)
   ```

3. **Linearized system:**
   ```
   λ_lin = -α - β·x_avg·⟨k⟩ + correction_term
   ```

---

## Data Flow

```
Network Data (MAT/TXT)
        ↓
Graph Construction (networkx)
        ↓
Adjacency Matrix (A), Degree Vector
        ↓
   ↙           ↓          ↘
[Analytical]  [Numerical]  [Visualization]
  Eigenvalue  ODE Solver    Phase Space
  Computation (solve_ivp)   Plots
        ↓          ↓          ↓
  [Results] ← Comparison → [Figures]
```

---

## Python Dependencies

- **networkx** - Network analysis and generation
- **numpy** - Numerical computing
- **scipy.integrate** - ODE solving (solve_ivp)
- **scipy.io** - MATLAB file I/O
- **scipy.linalg** - Eigenvalue computation (eig)
- **scipy.sparse.linalg** - Sparse eigenvalue computation (eigsh)
- **matplotlib.pyplot** - Visualization

---

## Key Research Questions

1. **How does network structure (topology) affect stability?**
   - Different networks: Erdős-Rényi, Barabási-Albert, etc.

2. **What is the relationship between eigenvalues and dynamical stability?**
   - Analytical vs. numerical eigenvalues

3. **How do heterogeneous degrees influence instability?**
   - Mean-field vs. degree-dependent predictions

4. **Where do bifurcations occur in parameter space?**
   - Critical (β, χ) values where λ_max = 0

---

## Interpretation of Results

### Bifurcation Points (Phase Transitions)

- **Below bifurcation:** System converges to stable endemic equilibrium
- **At bifurcation:** System undergoes instability transition
- **Above bifurcation:** System behavior changes fundamentally

### Network Effects

- **Scale-free networks (BA):** Different critical parameters than random (ER)
- **Degree heterogeneity:** Increases instability threshold
- **Network density:** Affects eigenvalue magnitude

---

## Usage Instructions

### Run Dynamics Simulation
```python
python SIS_FX_PT.py
```
Outputs plot of infection equilibrium vs. transmission rate β

### Analyze Stability
```python
python dynimics_fx_pt_stability.py
```
Generates eigenvalue comparison for analytical vs. numerical methods

### Create Visualizations
```python
python fig_1.py    # Phase space bifurcation plot
python fig_2(1).py # Eigenvector scatter plot
python fig_2(2).py # Network visualization
python fig_3.py    # Extended stability analysis
```

---

## Contact & Attribution

**Author:** Shraosi (based on file headers)

**Research Focus:** Dynamical systems stability on networks, epidemic modeling, spectral analysis
