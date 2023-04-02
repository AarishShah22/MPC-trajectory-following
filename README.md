# Model Predictive Control (MPC) for Trajectory Following

This repository contains the implementation of MPC for trajectory following using a kinematic bicycle model for a car.

## Model

We use the kinematic bicycle model for a car to implement MPC and simulate the vehicle.

$$ \begin{bmatrix} \dot{x} \\ \dot{y} \\ \dot{z} \end{bmatrix} = \begin{bmatrix} u cos(\psi)-\frac{b}{L}u tan(\delta) sin(\psi) \\ u sin(\psi)+\frac{b}{L}u tan(\delta) cos(\psi) \\ \frac{u}{L} tan(\delta) \end{bmatrix}  $$

where $(x,y,z)$ are the states (pose) of the vehicle and $(u,\delta)$ are the inputs (longitudinal velocity, wheel angle).

## MPC

We use the a time step of 0.01 seconds and a time horizon of 10 timesteps. We use the `quadprog` command in MATLAB for generating the control inputs.