---
title: Apply Transforms and Spectral Methods
description: Apply symbolic Fourier, Laplace, or Z-transforms to derive frequency-domain representations and convert to numeric objects for engineering analysis.
tags: [matlab, symbolic-math, laplace, ztrans, fourier, ilaplace, transforms, live-script]
release: All releases
notes:
---

# Apply Transforms and Spectral Methods

Uses Symbolic Math Toolbox to transform time-domain equations or signals into the frequency or z-domain, then applies operations like partial-fraction expansion and inverse transforms to obtain closed-form results. Useful for filter design, transfer-function derivation, pulse-shaping analysis, or discretizing continuous-time controllers for embedded implementation.

## Metadata

- **Tags:** `matlab` `symbolic-math` `laplace` `ztrans` `fourier` `ilaplace` `transforms` `live-script`
- **MATLAB Release:** All releases
- **Required Toolboxes:** Symbolic Math Toolbox, Control System Toolbox (if producing tf or ss objects)

## The Prompt

```text
Create a concise plain‑text MATLAB Live Code file (.m) showing symbolic
[fourier | laplace | ztrans] analysis of [SYSTEM / SIGNAL DESCRIPTION]
for [ENGINEERING GOAL].
--- [SIGNAL | SYSTEM | COMPENSATOR] DEFINITION ---
  [DEFINE THE TIME‑DOMAIN OR FREQUENCY‑DOMAIN EXPRESSIONS]
  [DEFINE ANY GOVERNING ODES IF STARTING FROM EQUATIONS OF MOTION]
--- PARAMETER VALUES ---
  [LIST EACH PARAMETER WITH VALUE AND UNITS]
Requirements
Use Symbolic Math Toolbox to:
1. Define the symbolic expression(s); apply [fourier() | laplace() | ztrans() or symbolic Tustin substitution] to obtain the transformed representation; simplify and display
2. Apply additional symbolic operations as appropriate:
   - partfrac() for partial‑fraction expansion
   - ifourier() / ilaplace() / iztrans() for inverse transform
   - collect() or coeffs() to extract polynomial coefficients
   Display results at each step
3. Substitute parameter values; convert to a numeric [tf | tf with Ts | FIR vector | coefficient arrays] using double() on the symbolic coefficients; display
4. Generate verification plots comparing symbolic and numeric results:
   [BODE PLOT | FFT OVERLAY | CONTINUOUS vs. DISCRETE BODE COMPARISON | TIME‑DOMAIN RESPONSE]
5. [OPTIONAL] Derive the time‑domain implementation form (difference equation, FIR filter coefficients, etc.) for embedded deployment; display
6. Save the symbolic expressions and numeric objects to a MAT‑file
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

### Example 1: Quarter-Car Suspension with Laplace

```
Create a concise plain‑text MATLAB Live Code file (.m) showing symbolic Laplace‑domain analysis of a 2‑DOF quarter‑car suspension model for ride comfort and handling evaluation.

--- MODEL DEFINITION ---
States:  x = [zs; dzs; zu; dzu]
         (sprung displacement, sprung velocity,
          unsprung displacement, unsprung velocity)
Input:   zr(t)  (road profile displacement)
Outputs: body acceleration ddzs (ride comfort),
         tire deflection zu - zr (handling / road holding)

Equations of motion:
  ms*ddzs = -ks*(zs - zu) - cs*(dzs - dzu)
  mu*ddzu =  ks*(zs - zu) + cs*(dzs - dzu) - kt*(zu - zr)
where ms is the sprung (body) mass, mu is the unsprung (wheel) mass, ks and cs are the suspension stiffness and damping, and kt is the tire stiffness.

--- PARAMETER VALUES ---
  ms  = 320    kg       (quarter‑body mass)
  mu  = 40     kg       (wheel + axle mass)
  ks  = 18000  N/m      (spring stiffness)
  cs  = 1500   N*s/m    (damper coefficient)
  kt  = 200000 N/m      (tire stiffness)

Requirements
Use Symbolic Math Toolbox to:
1. Apply laplace() to both ODEs (assuming zero initial conditions) and solve for Zs(s) and Zu(s) as functions of Zr(s).
   Then form the symbolic transmissibility for body acceleration per road displacement (H_body(s)), and for tire deflection deflection(H_tire(s))
   Simplify and display both
2. Apply the partial-fraction expansion to H_body(s) and display it
3. Compute the symbolic impulse response of H_body using ilaplace(), then display it and plot it using numeric parameters
4. Substitute parameter values
   Convert to numeric tf objects using double() on coefficients, and display both transfer functions
5. Generate a Bode magnitude plot of H_body and H_tire on the same axes (frequency in Hz).
   Annotate the body‑bounce and wheel‑hop resonant frequencies
6. Save the symbolic and numeric transfer functions to a MAT‑file
```

### Example 2: Discretize Compensator via Z-Transform

```
Create a concise plain‑text MATLAB Live Code file (.m) showing symbolic Tustin discretization of a Type II compensator for digital voltage‑loop control of a DC‑DC converter.

--- COMPENSATOR DEFINITION ---
Continuous‑time Type II compensator:
  C(s) = Kc * (1 + s/wz)^2 / (s * (1 + s/wp))
This provides two zeros for phase boost, one integrator for DC
regulation, and a high‑frequency pole for noise roll‑off.
Design targets (already met by the continuous design):
  Crossover frequency  fc  = 20 kHz
  Phase margin              ≈ 55°
Tustin (bilinear) transform:
  s  =  (2/Ts) * (z - 1) / (z + 1)

--- PARAMETER VALUES ---
  Kc   = 2.5e4           (compensator gain)
  fz   = 3.5e3  Hz       (zero frequency, wz = 2*pi*fz)
  fp   = 120e3  Hz       (pole frequency, wp = 2*pi*fp)
  fsw  = 200e3  Hz       (switching / sample frequency)
  Ts   = 1/fsw           (sample period)

Requirements
Use Symbolic Math Toolbox to:
1. Define C(s) symbolically, display it, and generate a Bode plot showing crossover and phase margin
2. Apply Tustin substitution symbolically to get C(z), then simplify it and display as a ratio of polynomials in z
3. Collect numerator and denominator coefficients of C(z) in descending powers of z, then display them as vectors b and a
4. Derive a time‑domain difference equation suitable for C implementation on a DSP and display it
5. Substitute parameter values.
   Convert to a numeric tf object with sample time Ts
   Overlay the discrete Bode plot against the continuous one to verify frequency‑domain match up to fsw/2
6. Save the numeric b/a coefficient vectors, the tf objects (continuous and discrete), and the sample rate to a MAT‑file
```

## Common Patterns

```matlab
% Pattern 1: Transforms
syms t s x w n z

% Laplace transform
laplace(exp(-2*t), t, s)             % 1/(s + 2)
laplace(sin(3*t), t, s)             % 3/(s^2 + 9)
laplace(t^2*exp(-t), t, s)          % 2/(s + 1)^3

% Inverse Laplace
ilaplace(1/(s^2 + 4), s, t)         % sin(2*t)/2

% Fourier transform
fourier(exp(-t^2), t, w)

% Inverse Fourier
ifourier(2*exp(-abs(w)), w, t)

% Z-transform
ztrans(2^n, n, z)                    % z/(z - 2)
iztrans(z/(z - 2), z, n)            % 2^n

% Pattern 2: Deriving Transfer Functions from Differential Equations
%% Transfer function of a mass-spring-damper system
% m*y'' + c*y' + k*y = F(t)
syms y(t) F(t) s Y U m c k

% Define the ODE
ode = m*diff(y,t,2) + c*diff(y,t) + k*y == F;

% Take Laplace transform
ode_laplace = laplace(ode, t, s);

% Substitute Laplace-domain variables and zero initial conditions
ode_s = subs(ode_laplace, ...
    [laplace(y(t), t, s), laplace(F(t), t, s), ...
     subs(diff(y(t), t), t, 0), y(0)], ...
    [Y, U, 0, 0]);

% Solve for Y in terms of U
Y_sol = solve(ode_s, Y);

% Transfer function G(s) = Y(s)/U(s)
G = simplify(Y_sol / U)
% Result: 1/(k + c*s + m*s^2)

%% Convert to numeric tf object (requires Control System Toolbox)
% Substitute numeric parameter values first
G_num = subs(G, [m, c, k], [1, 2, 5]);
[num_coeffs, den_coeffs] = numden(G_num);
num_poly = sym2poly(num_coeffs);
den_poly = sym2poly(den_coeffs);
tf_sys = tf(num_poly, den_poly);

% Pattern 2: Laplace Transform for Control Systems Analysis
%% System: RLC circuit
% L*di/dt + R*i + (1/C)*integral(i dt) = V(t)
% In terms of charge q: L*q'' + R*q' + q/C = V(t)

syms t s q(t) V(t) L R C Q Vs

% Define ODE
ode = L*diff(q,t,2) + R*diff(q,t) + q/C == V;

% Take Laplace transform
ode_L = laplace(ode, t, s);

% Substitute for zero ICs and Laplace variables
ode_s = subs(ode_L, ...
    [laplace(q(t),t,s), laplace(V(t),t,s), q(0), subs(diff(q(t),t),t,0)], ...
    [Q, Vs, 0, 0]);

% Solve for transfer function G(s) = Q(s)/Vs(s)
Q_sol = solve(ode_s, Q);
G = simplify(Q_sol / Vs)
% Result: 1/(1/C + R*s + L*s^2)

% Substitute numeric values and create tf
G_num = subs(G, [L, R, C], [1e-3, 100, 1e-6]);
[n, d] = numden(G_num);
sys = tf(sym2poly(n), sym2poly(d));
bode(sys)
```

## Related Prompts

- [Linearize Nonlinear Dynamics](./linearize.md)
- [Diagram to Symbolic Model](./diagram-to-model.md)
- [Symbolic Parameters in Numerical Workflows](./parametric.md)
- [Generate Plain Text Live Scripts](../live-scripts-documentation/generate-plain-text-live-script.md)

## References

- [MATLAB Documentation: Solve Differential Equations of RLC Circuit Using Laplace Transform](https://www.mathworks.com/help/symbolic/solve-differential-equations-using-laplace-transform.html)
- [MATLAB Documentation: Solve Difference Equations Using Z-Transform](https://www.mathworks.com/help/symbolic/compute-z-transforms-and-inverse-z-transforms.html)
- [MATLAB Documentation: Live Code File Format (.m)](https://www.mathworks.com/help/matlab/matlab_prog/plain-text-file-format-for-live-scripts.html)
- [Video: What are Transfer Functions?](https://www.youtube.com/watch?v=2Xl7--Df3g8)
- [Video: Understanding the Z-Transform](https://www.youtube.com/watch?v=XJRW6jamUHk)
---