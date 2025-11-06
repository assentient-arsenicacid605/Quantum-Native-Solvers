# Quantum-Native-Algos
A compilation of quantum-native solver techniques that can be used to map and run on a quantum computer. Compiled by Onri Jay Benally. 

The purpose of this repository is to monitor computational techniques, over time, that can be used to translate a working mathematical or computational framework of interest that may be eligible for quantum simulation. 

---



```
Quantum‑Native Simulation Techniques 
├─ I. Lattice‑kinetic/ Automata family
│   ├─ A. Quantum Lattice Boltzmann (QLBM)
│   │   ├─ A1. “QW‑equivalent” formulations (QLBM ≈ quantum walk)
│   │   ├─ A2. Carleman‑linearized LB (LB–Carleman, fault‑tolerant oriented)
│   │   ├─ A3. Multi‑circuit QLBM for noisy devices
│   │   └─ A4. Early hardware demos (advection–diffusion on QPU)
│   ├─ B. Quantum Lattice Gas Automata (QLGA)
│   │   ├─ B1. Measurement‑based QLG for Navier–Stokes
│   │   ├─ B2. Efficient QLGA single‑step schemes
│   │   └─ B3. Fully‑quantum LGA building blocks (MHD/NS‑ready)
│   └─ C. Quantum Walks/ Quantum Cellular Automata (QCA)
│       ├─ C1. Dirac‑type dynamics (transport, diffusion, barriers & tunneling)
│       ├─ C2. Schrödinger‑type dynamics (continuous‑time QW ↔ free‑particle propagation,
│       │      discrete‑time QW approximating Laplacian)
│       └─ C3. Gauge‑field‑coupled QW/QCA (U(1), SU(2) links for both Dirac & Schrödinger)
│
├─ II. Quantum PDE & Linear‑Differential Equation solvers
│   ├─ D. Quantum Linear‑System Algorithms (HHL & block‑encoding, Poisson)
│   ├─ E. Linear Differential Equations (LDE) oracles (Berry–Childs–Kothari)
│   ├─ F. Carleman linearization for nonlinear ODE/PDE (Navier–Stokes regimes)
│   ├─ G. Quantum Finite Elements (Qu‑FEM/ Q‑FEM) + q‑Multigrid
│   ├─ H. Dirac‑PDE solvers
│   │   ├─ H1. QW‑based Dirac first‑order hyperbolic solvers
│   │   ├─ H2. QLGA spinor kernels for relativistic fluids
│   │   └─ H3. LCU/ qubitization of the Dirac operator (mass + kinetic terms)
│   └─ I. Schrödinger‑PDE solvers
│       ├─ I1. Spectral (Fourier) discretization + Trotter/LCU propagation
│       ├─ I2. Chebyshev‑polynomial and interaction‑picture schemes
│       └─ I3. Variational linear‑solver (VQLS) for stationary Schrödinger equation
│
├─ III. Variational & Hybrid Dynamics
│   ├─ J. Variational Quantum Simulation (VQS: real & imaginary time)
│   ├─ K. Quantum Imaginary‑Time Evolution (QITE) + accelerators
│   ├─ L. Variational Quantum Linear Solver (VQLS) for scattering/PDE
│   ├─ M. Quantum Car‑Parrinello Molecular Dynamics (QCPMD, NISQ)
│   ├─ N. Variational Dirac solvers
│   │   ├─ N1. Real‑time VQS for Dirac spinor dynamics
│   │   └─ N2. QITE for massive Dirac ground‑state preparation
│   └─ O. Variational Schrödinger solvers
│       ├─ O1. Real‑time VQS with hardware‑efficient ansätze for non‑relativistic dynamics
│       └─ O2. Imaginary‑time VQS/ QITE for ground‑state Schrödinger problems
│
├─ IV. Hamiltonian Simulation & Scattering
│   ├─ P. Real‑space electron‑dynamics (LCU/ qubitization/ Trotter)
│   ├─ Q. Scattering & cross‑section algorithms (Green’s functions, S‑matrix)
│   ├─ R. Quantum‑walk tunneling/ barrier transport
│   ├─ S. Dirac‑Hamiltonian simulation
│   │   ├─ S1. Product‑formula (Trotter‑Suzuki) for γ·p + m terms
│   │   ├─ S2. LCU/ qubitization of relativistic kinetic operator
│   │   └─ S3. Interaction‑picture and Dyson‑series methods for time‑dependent Dirac fields
│   └─ T. Schrödinger‑Hamiltonian simulation
│       ├─ T1. Finite‑difference/plane‑wave discretization + Trotter
│       ├─ T2. LCU/ qubitization of kinetic + potential operators
│       └─ T3. Interaction‑picture schemes for driven potentials
│
├─ V. Stochastic & Monte Carlo family
│   ├─ O. Quantum Monte Carlo (classically executed, physics‑native)
│   │   ├─ O1. Variational Monte Carlo (VMC; |Ψ|² sampling)
│   │   ├─ O2. Diffusion Monte Carlo (DMC; imaginary‑time projection, fixed‑node)
│   │   ├─ O3. Path‑Integral QMC (PIMC; finite‑T & ground‑state PIGS)
│   │   ├─ O4. Auxiliary‑Field QMC (AFQMC; Hubbard–Stratonovich)
│   │   └─ O5. Stochastic Series Expansion (SSE; loop/cluster updates)
│   ├─ P. Quantum‑accelerated Monte Carlo (hardware‑native)
│   │   ├─ P1. Quantum Amplitude Estimation (QAE; quadratic sampling speedups)
│   │   ├─ P2. Quantum Metropolis sampling (thermal states without sign problem)
│   │   └─ P3. Hybrid QMC w/ QC‑generated trials (e.g., QSCI → AFQMC)
│   └─ Q. Quantum Trajectories (Monte Carlo wave function, stochastic Schrödinger)
│       ├─ Q1. Quantum jumps/ MCWF unraveling of Lindblad
│       └─ Q2. Adiabatic‑master‑equation trajectories; restricted‑MCWF diagnostics
│   ├─ R. Schrödinger‑QMC variants
│   │   ├─ R1. Diffusion Monte Carlo with quantum‑enhanced trial‑wave‑function optimization
│   │   └─ R2. Path‑Integral QMC for real‑time (Keldysh) Schrödinger evolution
│   └─ S. Dirac‑QMC variants
│       ├─ S1. Stochastic Dirac equation (stochastic quantization of relativistic fields)
│       └─ S2. Quantum‑walk based Monte Carlo for relativistic barrier tunneling
│
├─ VI. Open‑System (Lindblad) solvers & Magnetization links
│   ├─ R. Lindblad simulators on quantum hardware
│   │   ├─ R1. LCU/ Stinespring‑dilation & qubitization‑style algorithms
│   │   ├─ R2. Randomized product‑formula & sampling compilers for GKSL
│   │   └─ R3. Quantum‑trajectory simulators with additive O(T + log 1/ε) scaling (class‑dependent)
│   ├─ S. Classical‑HPC Lindblad solvers with quantum trajectories (large‑N pure‑state sampling)
│   └─ T. Lindblad → LL/LLG mappings for magnetization dynamics
│       ├─ T1. LL from Lindbladian dynamics (scale‑separated local mean‑field regimes)
│       ├─ T2. Microscopic LLG coefficients via NEGF/scattering formulations
│       └─ T3. Quantum‑LLG analogs and constraints from entanglement/ non‑Markovianity
│
├─ VII. Emerging quantum‑native simulation families
│   ├─ A. Quantum Tensor‑Network Simulators (MPS, PEPS, TTN) for PDEs & many‑body states
│   ├─ B. Quantum Signal Processing/ QSVT‑based PDE solvers
│   │        (Chebyshev and polynomial transformations for evolution and linear systems)
│   ├─ C. Fractional & Non‑local Operator Solvers
│   │        (fractional Laplacian, Riesz derivatives, etc.)
│   ├─ D. Quantum Stochastic Differential Equations (QSDE) &
│   │        non‑Markovian open‑system dynamics
│   ├─ E. Quantum Lattice‑Gauge‑Theory Engines
│   │        (U(1), SU(2), SU(3) gauge fields, plaquette unitaries)
│   ├─ F. Quantum‑Accelerated Adaptive‑Mesh‑Refinement (AMR)
│   │        (amplitude‑estimation‑driven error indicators)
│   ├─ G. Quantum‑Enhanced Real‑Time Time‑Dependent Density‑Functional Theory (RT‑TDDFT)
│   ├─ H. Quantum Wigner‑Function/ Phase‑Space Solvers
│   │        (discrete Wigner representation + FFT‑based evolution)
│   ├─ I. Quantum Uncertainty Quantification (UQ) for PDEs
│   │        (parameter‑sampling via QAE)
│   ├─ J. Quantum Spectral‑Element/ hp‑Finite‑Element Methods
│   │        (high‑order basis encoded in qubit registers)
│   ├─ K. Quantum Metropolis‑Hastings for Thermodynamic Ensembles of PDE Solutions
│   ├─ L. Hybrid Quantum‑Classical Domain Decomposition (Schwarz) for large‑scale PDEs
│   ├─ M. Quantum‑Accelerated Sparse‑Matrix Preconditioners
│   │        (quantum BiCGSTAB, quantum multigrid)
│   └─ N. Quantum Embedding & Multi‑Scale Methods
│           (DMFT, QM/MM, impurity solvers, quantum‑accelerated self‑energy evaluation)
│
├─ VIII. Dirac Equation Solvers (relativistic quantum dynamics)
│   ├─ A. Quantum‑Walk‑based Dirac solvers
│   │   ├─ A1. 1‑D split‑step Dirac QW (continuous limit → Dirac equation)
│   │   ├─ A2. 2‑D/ 3‑D multi‑step Dirac QWs (including spin‑orbit coupling)
│   │   ├─ A3. Gauge‑field‑coupled Dirac QWs (U(1) electromagnetic, SU(2) weak)
│   │   └─ A4. Multi‑circuit Dirac QWs for NISQ devices (error‑mitigated implementations)
│   ├─ B. Quantum Lattice‑Gas Automata for Dirac dynamics
│   │   ├─ B1. Spinor QLGA kernels (γ‑matrix encoding on qubits)
│   │   ├─ B2. Dirac–Maxwell QLGA (self‑consistent coupling of spinor and EM fields)
│   │   └─ B3. Relativistic fluid‑LGA (Dirac‑MHD, relativistic Navier–Stokes)
│   ├─ C. Hamiltonian‑simulation of the Dirac operator
│   │   ├─ C1. Trotter‑Suzuki product formulas for γ·p + m terms
│   │   ├─ C2. LCU/ qubitization of the Dirac Hamiltonian (including mass‑gap)
│   │   ├─ C3. Interaction‑picture/ Dyson‑series Dirac simulators for time‑dependent fields
│   │   └─ C4. Spectral‑method (Fourier) Dirac simulation (momentum‑space discretization)
│   ├─ D. Variational & Imaginary‑time Dirac solvers
│   │   ├─ D1. Real‑time VQS for Dirac spinors (hardware‑efficient ansätze)
│   │   ├─ D2. QITE for massive Dirac ground‑state preparation
│   │   └─ D3. Adaptive‑VQE style Dirac dynamics (ADAPT‑VQD)
│   ├─ E. Carleman‑linearized Dirac‑type nonlinear models
│   │   ├─ E1. Dirac–Klein‑Gordon coupling (scalar‑field interactions)
│   │   ├─ E2. Non‑linear Dirac (Thirring, Gross‑Neveu) via Carleman embedding
│   │   └─ E3. Dirac–Maxwell–Higgs hybrid systems
│   ├─ F. Quantum Finite‑Element/ Multigrid for Dirac
│   │   ├─ F1. QFEM‑Dirac (mesh‑based spinor discretization)
│   │   └─ F2. Dirac‑multigrid preconditioners realized with quantum subroutines
│   ├─ G. Scattering & S‑matrix for relativistic particles
│   │   ├─ G1. Green‑function based Dirac cross‑section algorithms
│   │   ├─ G2. Relativistic quantum‑walk tunneling analysis (Klein paradox)
│   │   └─ G3. Dirac‑time‑dependent scattering (time‑varying potentials)
│   └─ H. Hybrid classical‑quantum Dirac workflows
│       ├─ H1. Classical spectral Dirac solver + quantum subroutines for high‑frequency modes
│       └─ H2. Domain‑decomposition: QC solves interior relativistic region, classical for exterior
│
└─ IX. Schrödinger Equation Solvers (non‑relativistic quantum dynamics)
    ├─ A. Direct Hamiltonian‑simulation of Schrödinger dynamics
    │   ├─ A1. Spatial discretization (finite‑difference, plane‑wave) + Trotter‑Suzuki
    │   ├─ A2. LCU/ qubitization of kinetic (∇²) + potential V(x) operators
    │   ├─ A3. Interaction‑picture/ Dyson‑series approaches for time‑dependent potentials
    │   └─ A4. High‑order product formulas (e.g., Forest‑Ruth, Suzuki‑23) for improved accuracy
    ├─ B. Quantum‑PDE solvers specialized for Schrödinger
    │   ├─ B1. Linear‑Differential‑Equation oracle tailored to time‑dependent Schrödinger equation
    │   ├─ B2. Spectral (Fourier) method via Quantum Fourier Transform (QFT) – plane‑wave basis
    │   ├─ B3. Chebyshev‑polynomial expansion of the propagator e^{-iHt}
    │   └─ B4. Pseudo‑spectral split‑operator schemes (FFT‑based) on quantum hardware
    ├─ C. Variational & Hybrid Schrödinger solvers
    │   ├─ C1. Real‑time VQS with hardware‑efficient ansätze (e.g., hardware‑efficient, ADAPT‑VQE)
    │   ├─ C2. Imaginary‑time VQS/ QITE for ground‑state preparation of Schrödinger Hamiltonians
    │   ├─ C3. Variational Quantum Linear Solver (VQLS) for stationary Schrödinger equation (Ax = b with A = H−E)
    │   └─ C4. Adaptive‑VQD for excited‑state dynamics
    ├─ D. Quantum‑Walk implementations of Schrödinger dynamics
    │   ├─ D1. Continuous‑time quantum walk mapping to free‑particle Schrödinger (H = −Δ)
    │   ├─ D2. Discrete‑time quantum walk reproducing Laplacian kinetic term via coin‑shift construction
    │   └─ D3. Quantum‑walk based tunneling simulations for arbitrary barrier profiles
    ├─ E. Quantum Cellular Automata for Schrödinger equation
    │   ├─ E1. Local unitary update rules approximating the Laplacian operator
    │   └─ E2. Multi‑step QCA with on‑site potential encoding
    ├─ F. Quantum Monte Carlo methods for Schrödinger equation
    │   ├─ F1. Diffusion Monte Carlo (imaginary‑time projection) with QC‑accelerated trial‑wave‑function refinement
    │   ├─ F2. Path‑Integral QMC for real‑time evolution (Keldysh contour) using quantum‑enhanced sampling
    │   ├─ F3. Variational Monte Carlo with quantum‑generated ansätze (QANSAT)
    │   └─ F4. Reptation QMC with quantum‑assisted move proposals
    ├─ G. Many‑body Schrödinger solvers (second‑quantized)
    │   ├─ G1. Electronic‑structure Hamiltonian simulation (quantum chemistry) via LCU/ qubitization
    │   ├─ G2. Quantum Dynamical Mean‑Field Theory (QDMFT) for lattice models
    │   ├─ G3. Real‑time Green‑function methods (Kadanoff‑Baym) on quantum processors
    │   └─ G4. Tensor‑network‑inspired hybrid algorithms (e.g., VQE‑DMRG hybrids)
    ├─ H. Hybrid Classical‑Quantum Schrödinger solvers
    │   ├─ H1. Classical spectral (Fourier) solver + quantum subroutine for high‑frequency phase propagation
    │   ├─ H2. Domain‑decomposition: QC handles sub‑grid potentials; classical handles bulk
    │   └─ H3. Multi‑scale coupling of quantum‑accelerated short‑range dynamics with classical long‑range Poisson solver
    ├─ I. Scattering & S‑matrix for non‑relativistic particles
    │   ├─ I1. Green‑function based cross‑section computation (Lippmann‑Schwinger equation)
    │   ├─ I2. Quantum‑walk tunneling analysis for arbitrary 1‑D/ 2‑D barrier profiles
    │   └─ I3. Real‑time S‑matrix extraction via phase‑estimation on scattered wave‑packets
    └─ J. Quantum Finite‑Element/ Multigrid for Schrödinger
        ├─ J1. Mesh‑based Schrödinger solvers with quantum enforcement of continuity/ boundary conditions
        ├─ J2. Quantum‑accelerated multigrid preconditioners for the Laplacian‑plus‑potential operator
        └─ J3. Adaptive mesh refinement driven by quantum‑measured error estimates
```
