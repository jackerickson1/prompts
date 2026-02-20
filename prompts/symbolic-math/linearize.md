---
title: Linearize Nonlinear Dynamics for Control Design
description: Symbolically linearize nonlinear equations of motion around an operating point and export a numeric state-space model for control design
tags: [matlab, symbolic-math, control-design, linearization, jacobian, state-space, trim, live-script]
release: All releases / Requires Symbolic Math Toolbox, Control System Toolbox
notes:
---

# Linearize Nonlinear Dynamics for Control Design

Takes a set of nonlinear ODEs, uses Symbolic Math Toolbox to compute Jacobians and trim conditions, then substitutes parameter values to produce a linear state-space model compatible with Control System Toolbox. Useful whenever you have first-principles dynamics and need a linear plant model for control design.

## Metadata

- **Tags:** `matlab` `symbolic-math` `control-design` `linearization` `jacobian` `state-space` `trim` `live-script`
- **MATLAB Release:** All releases
- **Required Toolboxes:** Symbolic Math Toolbox, Control System Toolbox

## The Prompt

```text
Create a concise plain‑text MATLAB Live Code file (.m) showing
symbolic‑to‑numeric linearization of [SYSTEM DESCRIPTION] for
[CONTROL DESIGN GOAL].
--- MODEL DEFINITION ---
States:  x = [STATE VARIABLES AND DESCRIPTIONS]
Inputs:  u = [INPUT VARIABLES AND DESCRIPTIONS]
Nonlinear equations of motion:
  [dx1/dt = f1(x, u)  ]
  [dx2/dt = f2(x, u)  ]
  [  ...               ]
Auxiliary expressions:
  [FORCES, MOMENTS, CONSTITUTIVE RELATIONS, ETC.]
--- PARAMETER VALUES ---
  [LIST EACH PARAMETER WITH VALUE AND UNITS]
--- OPERATING POINT ---
  [SPECIFY THE EQUILIBRIUM CONDITION: either prescribe the state
   values directly, or state the constraint (e.g., "steady level
   flight: all derivatives = 0") and let solve find the trim point]

Requirements
Use Symbolic Math Toolbox to:
1. Define the nonlinear state equations symbolically; display them
2. Determine the operating point: evaluate or solve for equilibrium states and inputs; display the trim values
3. Compute symbolic Jacobians A = df/dx and B = df/du; display both
4. Substitute parameter values and operating point; convert to numeric with double(); display A and B and compute eigenvalues
5. Build a state‑space model ss(A,B,C,D) with appropriate C and D for [MEASURED OUTPUTS]; generate a Bode plot from [INPUT] to [OUTPUT]
6. Save the ss object and relevant symbolic expressions to a MAT‑file

[LIVE SCRIPT FORMATTING]
```

## Usage Tips

1. **Live Script Formatting:** Symbolic equations display best as Live Script output. Starting with R2025a, MATLAB supports plain text Live Code formatting. For best results, use the [MATLAB Plain Text Live Script Generator](https://github.com/matlab/skills/blob/main/skills/matlab-live-script/SKILL.md) skill. Alternatively, you can include these instructions:
```text
* Use disp() for all outputs (no fprintf, no pretty())
* Section heads should look like this:
    %%
    %[text] ## Section Title
* DO NOT do this: `%% Section Title`
* Keep comments brief
* Display: symbolic equations, trim point, symbolic and numeric A/B, eigenvalues, Bode plot
Output: Single executable script, minimal but complete, include the appendix metadata block at the end (mandatory for plain‑text Live Script recognition):
    %[appendix]{"version":"1.0"}
    %---
    %[metadata:view]
    %   data: {"layout":"inline"}
    %---
```
2. **Symbolic Math Toolbox skill:** For best results, load the [Symbolic Math Toolbox skill](https://github.com/matlab/skills/blob/main/skills/symbolic-math-toolbox/SKILL.md)
3. **Display equations and generate plots:** For each step, generate outputs so you can verify the generated code.

## Example Usage

### Example 1: Trim & Linearize for Aircraft Flight Control

```
Create a concise plain‑text MATLAB Live Code file (.m) showing symbolic-to-numeric linearization of aircraft longitudinal dynamics.

--- MODEL DEFINITION ---

States:  x = [V; alpha; q; theta]
         (airspeed, angle of attack, pitch rate, pitch angle)
Inputs:  u = [delta_e; delta_T]
         (elevator deflection, thrust)

Nonlinear equations of motion:

  dV/dt     = (T*cos(alpha) - D - m*g*sin(theta - alpha)) / m
  dalpha/dt = q - (T*sin(alpha) + L - m*g*cos(theta - alpha)) / (m*V)
  dq/dt     = M_pitch / Iy
  dtheta/dt = q

Aerodynamic forces and moment:

  qbar = 0.5 * rho * V^2
  L    = qbar * S * (CL0 + CLa*alpha + CLde*delta_e)
  D    = qbar * S * (CD0 + CDa*alpha)
  T    = delta_T
  M_pitch = qbar * S * cbar * (Cm0 + Cma*alpha + Cmq*(cbar/(2*V))*q + Cmde*delta_e)

--- PARAMETER VALUES ---

  m     = 1200   kg        (aircraft mass)
  Iy    = 2100   kg*m^2    (pitch inertia)
  S     = 18     m^2       (wing area)
  cbar  = 1.5    m         (mean aerodynamic chord)
  rho   = 1.225  kg/m^3    (air density, sea level)
  g     = 9.81   m/s^2
  CL0   = 0.28             CLa  = 4.8   CLde = 0.36
  CD0   = 0.03             CDa  = 0.35
  Cm0   = 0.05             Cma  = -0.5  Cmq = -8.0  Cmde = -1.1

Requirements
Use Symbolic Math Toolbox to:
1. Define symbolic equations of motion and display them symbolically
2. Trim for steady level flight and display the trim point
3. Compute symbolic Jacobians A and B, and display them
4. Substitute the parameter values and build a state-space model
5. Evaluate numerically and compute eigenvalues to show short-period and phugoid modes. Generate a Bode plot from elevator to pitch rate
6. Generate a MATLAB state‑space object ss(A,B,C,D) stored in a MAT‑file
```

### Example 2: Linearize a 2-Link Robot Arm for Joint-Space PD Control

```
Create a concise plain‑text MATLAB Live Code file (.m) showing symbolic Lagrangian derivation and linearization of a planar 2‑link robot arm for joint‑space control design.

--- MODEL DEFINITION ---

Generalized coordinates:  q = [q1; q2]   (joint angles)
States:  x = [q1; q2; dq1; dq2]
Inputs:  u = [tau1; tau2]                 (joint torques)

Link geometry (both links modeled as uniform rods, mass at center):
  Position of center of mass of link 1:
    x1 = Lc1*cos(q1)
    y1 = Lc1*sin(q1)
  Position of center of mass of link 2:
    x2 = L1*cos(q1) + Lc2*cos(q1+q2)
    y2 = L1*sin(q1) + Lc2*sin(q1+q2)
Kinetic energy:
  T = 0.5*m1*(dx1^2+dy1^2) + 0.5*I1*dq1^2
    + 0.5*m2*(dx2^2+dy2^2) + 0.5*I2*(dq1+dq2)^2
Potential energy:
  V = m1*g*y1 + m2*g*y2
Euler‑Lagrange equations:
  M(q)*ddq + C(q,dq)*dq + G(q) = tau
  where M is the mass matrix, C the Coriolis/centrifugal matrix,
  G the gravity vector. Derive these from T and V symbolically.

--- PARAMETER VALUES ---
  m1   = 3.0   kg       (link 1 mass)
  m2   = 2.0   kg       (link 2 mass)
  L1   = 0.4   m        (link 1 length)
  L2   = 0.3   m        (link 2 length)
  Lc1  = 0.2   m        (link 1 center of mass)
  Lc2  = 0.15  m        (link 2 center of mass)
  I1   = 0.04  kg*m^2   (link 1 inertia about center)
  I2   = 0.015 kg*m^2   (link 2 inertia about center)
  g    = 9.81  m/s^2

--- OPERATING POINT ---
  q1_0  = 0      rad    (link 1 horizontal)
  q2_0  = pi/2   rad    (link 2 vertical, elbow up)
  dq1_0 = 0      rad/s
  dq2_0 = 0      rad/s

Requirements

Use Symbolic Math Toolbox to:
1. Define symbolic Lagrangian dynamics for a 2-link planar arm, then derive and display M(q), C(q,dq), G(q) 
2. Evaluate G at the operating point to find the equilibrium torques that hold the arm stationary; display the trim torques 
3. Compute symbolic Jacobians A = df/dx and B = df/du; display both matrices 
4. Substitute parameter values and operating point into A and B; convert to numeric with double(); display both matrices and compute eigenvalues 
5. Build a state-space model with C = eye(4) (full state feedback available) then generate a Bode plot from tau1 to q1 
6. Save the ss object and the symbolic M, C, G expressions to a MAT file
```

## Common Patterns

```matlab
% Pattern 1: Eigenvalues and eigenvectors
syms a b c d
A = [a b; c d];
[V, D] = eig(A)

% Pattern 2: Jacobian
syms x y
f = [x^2*y; 5*x + sin(y)];
J = jacobian(f, [x, y])
% [2*x*y, x^2; 5, cos(y)]

% Pattern3: Variable-precision arithmetic (VPA)
% Default: 32 significant digits
vpa(pi)                          % 3.1415926535897932384626433832795

% Change precision
digits(50);
vpa(pi)                          % 50-digit pi

% IMPORTANT: Convert to symbolic FIRST, then use vpa
% WRONG: vpa(1/3) gives vpa of the double 0.333... (already rounded)
% RIGHT: vpa(sym(1)/sym(3)) gives exact 1/3 to desired digits

% Reset precision
digits(32);                      % Back to default
```

## Related Prompts

- [Symbolic Transforms](./transforms.md)
- [Diagram to Symbolic Model](./diagram-to-model.md)
- [Symbolic Parameters in Numerical Workflows](./parametric.md)
- [Design PID Controller](../control-systems/design-pid-controller.md)
- [Generate Plain Text Live Scripts](../live-scripts-documentation/generate-plain-text-live-script.md)

## References

- [MATLAB Documentation: Jacobian matrix](https://www.mathworks.com/help/symbolic/sym.jacobian.html)
- [MATLAB Documentation: Choose Numeric or Symbolic Arithmetic](https://www.mathworks.com/help/symbolic/choose-symbolic-or-numeric-arithmetic.html)
- [MATLAB Documentation: Live Code File Format (.m)](https://www.mathworks.com/help/matlab/matlab_prog/plain-text-file-format-for-live-scripts.html)
- [Video: What is Linearization?](https://www.youtube.com/watch?v=5gEattuH3tI)
---