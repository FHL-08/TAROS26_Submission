# Constrained Torque-Level Formation Control of Nonholonomic Wheeled Mobile Robots via Consensus-Based Feedback Linearised MPC

**Abstract:**     This paper presents a distributed, torque-level formation control framework for differential-drive wheeled mobile robots. While nonlinear MPC (NMPC) and feedback-linearised MPC methods can achieve accurate formation tracking, their torque-level implementations are hindered by non-convex constraints, which lead to high computational cost. We address this challenge by combining feedback linearisation with a distributed MPC architecture and introducing a tractable Second-Order Cone Programming (SOCP) reformulation of the resulting non-convex problem. The proposed approach captures actuator-level dynamics while convexifying state-dependent constraints, enabling efficient optimisation with formal guarantees of recursive feasibility and local asymptotic stability. We benchmark the controller against NMPC, Linear MPC, and non-convex feedback-linearised MPC. Simulations show that the SOCP achieves torque smoothness comparable to the non-convex formulation while reducing average solver time from tens of milliseconds to around 2 ms in MATLAB, making it suitable for real-time implementation on computationally constrained platforms. These results position the method as a practical and scalable foundation for torque-level, distributed formation control on constrained systems.

---

## Citation

If you use this code in your research, please cite the following paper:

```bibtex
@inproceedings{Lawan2026TAROS,
  author    = {Faisal Lawan, Joaquin Carrasco Gomez, Weishu Zhan, Wei Pan},
  title     = {Constrained Torque-Level Formation Control of Nonholonomic Wheeled Mobile Robots via Consensus-Based Feedback Linearised MPC},
  booktitle = {Proceedings of Towards Autonomous Robotic Systems (TAROS)},
  year      = {2026},
}
```

## Overview
This repository contains the MATLAB simulation code for the results presented in the TAROS 2026 submission titled "Constrained Torque-Level Formation Control of Nonholonomic Wheeled Mobile Robots via Consensus-Based Feedback Linearised MPC". The code implements and benchmarks four different distributed model predictive control strategies for the formation control of nonholonomic wheeled mobile robots operating at the torque level.

The controllers included are:
1. **LMPC**: Jacobian-based Linearised MPC.
2. **NMPC**: Nonlinear MPC using the full robot dynamics.
3. **NFLMPC**: Non-convex Feedback Linearised MPC.
4. **SOCP**: Feedback Linearised MPC reformulated as a Second-Order Cone Programme.

## Repository Structure
The repository is organised into folders corresponding to the different simulation scenarios discussed in the paper:
```
.
├── Scenario_1/
│   ├── APF_LMPC_version.m
│   ├── APF_NMPC_version.m
│   ├── nonconvex_optimised_version_APF.m
│   └── SOCP_optimised_version_APF.m
├── Scenario_2/
│   └── ... (similar files)
├── Scenario_3/
│   └── ... (similar files)
├── Scenario_4/
│   └── ... (similar files)
├── Scenario_5/
│   └── ... (similar files)
├── LICENSE
└── README.md
```
Each Scenario_X folder contains the MATLAB scripts to run the simulations for the four controllers under that specific scenario's conditions (e.g., obstacle positions).

## Dependencies
To run these simulations, you will need:
1. **MATLAB:** (Tested on R2023b).
2. **Required MATLAB Toolboxes:**
    - Control System Toolbox
    - Optimization Toolbox
3. **External MATLAB Toolboxes/Libraries:**
    - **CasADi:** Required for symbolic representation, automatic differentiation, and interfacing with NLP/QP solvers. (https://web.casadi.org/)
    - **YALMIP:** Required for modelling and solving the SOCP formulation. (https://yalmip.github.io/)
4. **Solvers:**
    - **IPOPT:** An NLP solver required for the NMPC and NFLMPC controllers. (included with CasADi)
    - **MOSEK:** A high-performance solver for convex optimisation, including SOCP. Recommended for the SOCP controller. (https://www.mosek.com/)
    - **qpOASES:** A QP solver required for the LMPC controller. (included with CasADi)

Please ensure all required toolboxes are installed and added to the MATLAB path.

## Usage
1. Navigate to the desired `Scenario_X` folder within MATLAB.
2. Run the MATLAB script corresponding to the controller you wish to simulate (e.g.`SOCP_optimised_version_APF.m`).
3. The script will execute the simulation, display the key performance metrics in the command window, and generate a series of plots showing the results.
4. Simulation data is saved to `APF_Trajectory_Values.mat` in the respective scenario folder upon completion for further analysis.

## License
This project is licensed under the MIT License. See the [LICENSE](LICENSE) file for details.

Copyright &copy; 2026 Faisal Lawan.

## Contact
For questions regarding the paper or the code, please refer to the contact information provided in the conference paper.
