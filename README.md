# RF-MPC

Representation-Free Model Predictive Control (RF-MPC) is a MATLAB simulation framework for dynamic legged robots. RF-MPC represents the orientation using the rotation matrix and thus does not have the singularity issue associated with the Euler angles. The linear dynamics on the rotation matrix is derived using variation-based linearization (VBL).

## Requirement
Basic: MATLAB and MATLAB optimization toolbox
Optional: qpSWIFT (an efficient QP solver to MPC problems)

## Installation
There is no need to install external packages.

## Usage
navigate to the root directory and run the MAIN.m function

'''MATLAB
MAIN
'''

## Contributing
For major changes, please open an issue first to discuss what you would like to change.
