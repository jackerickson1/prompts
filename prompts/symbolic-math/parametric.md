---
title: Parameterize Models for Numerical Sweeps and Optimization
description: Keep design or operating parameters symbolic through derivation, then export callable MATLAB functions for parameter sweeps, optimization, or gain scheduling.
tags: [matlab, symbolic-math, parametric, optimization, matlabFunction, live-script]
release: All releases / Requires Symbolic Math Toolbox
notes:
---

# Parameterize Models for Numerical Sweeps and Optimization

Uses Symbolic Math Toolbox to derive closed-form expressions where selected parameters remain symbolic while others are substituted numerically, then exports optimized function files via `matlabFunction` with the 'File' option. Useful when an engineer needs to rapidly evaluate hundreds of design candidates or map performance across an operating envelope without re-deriving the model each time.

## Metadata

- **Tags:** `matlab` `symbolic-math` `parametric` `optimization` `matlabFunction` `live-script`
- **MATLAB Release:** All releases
- **Required Toolboxes:** Symbolic Math Toolbox, plus appropriate application toolboxes

## The Prompt

```text
Create a concise plain‑text MATLAB Live Code file (.m). The [FIXED CONTEXT: e.g., "vehicle spec"] sets [FIXED PARAMETERS], but [DESIGN / OPERATING VARIABLES] are being selected or vary across operating conditions. Export callable MATLAB functions that accept ([SYMBOLIC VARIABLE LIST]) and return [PERFORMANCE METRICS / SYSTEM OBJECTS], for use in [DOWNSTREAM WORKFLOW: optimization loop, Simulink parameter sweep, gain‑scheduling table, risk dashboard, etc.].

--- MODEL DEFINITION ---
  [GOVERNING EQUATIONS, TRANSFER FUNCTIONS, OR ANALYTIC FORMULA]

--- FIXED PARAMETERS (substitute numeric) ---
  [LIST EACH WITH VALUE AND UNITS]

--- DESIGN / OPERATING VARIABLES (keep symbolic) ---
  [LIST EACH WITH DESCRIPTION AND SWEEP RANGE]

Requirements
Use Symbolic Math Toolbox to:
1. Derive the symbolic model with [DESIGN VARIABLES] as syms and [FIXED PARAMETERS] already substituted numerically; simplify and display; show how key expressions depend on the symbolic variables.
2. Derive symbolic performance metrics or characteristic quantities (e.g., natural frequencies, damping ratios, sensitivities, Greeks); display.
3. Use matlabFunction() with the 'File' option to export one or more .m function files to disk that accept ([SYMBOLIC VARIABLE VALUES]) and return [NUMERIC OUTPUTS: tf objects, metric values, matrices, etc.]; verify at a representative test point against a manually-constructed result.
4. Sweep the design / operating variables over the specified ranges using the exported functions; compute the performance metrics at each grid point.
5. Generate design‑space visualizations:
   [CONTOUR / SURFACE PLOTS OF METRICS vs. VARIABLES | PARETO FRONT SHOWING TRADE‑OFFS | BOUNDARY OVERLAY SHOWING REQUIREMENT COMPLIANCE]
6. Save symbolic expressions, sweep data, and visualization‑ready results to a MAT‑file.
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

### Example 1: Aerodynamic Linearization and Gain Scheduling

```
Create a concise plain‑text MATLAB Live Code file (.m) that keeps key aerodynamic parameters symbolic through linearization and stability analysis, then exports numeric functions for gain scheduling across a flight envelope.

--- MODEL DEFINITION ---
Simplified short‑period approximation (2 states):
  States:  x = [alpha; q]   (angle of attack, pitch rate)
  Input:   u = delta_e       (elevator)
  State equations:
    dalpha/dt = Zalpha*alpha + q + Zde*delta_e
    dq/dt     = Malpha*alpha + Mq*q + Mde*delta_e
  where Zalpha, Zde, Malpha, Mq, Mde are dimensional stability and   control derivatives that vary with flight condition.
The A and B matrices are:
  A = [Zalpha, 1; Malpha, Mq]
  B = [Zde; Mde]

--- PARAMETER RELATIONSHIPS ---
Dimensional derivatives depend on dynamic pressure qbar and airspeed V:
  Zalpha = qbar*S*CLa / (m*V)      (negative, restoring)
  Malpha = qbar*S*cbar*Cma / Iy
  Mq     = qbar*S*cbar^2*Cmq / (2*Iy*V)
  Zde    = qbar*S*CLde / (m*V)
  Mde    = qbar*S*cbar*Cmde / Iy
  qbar   = 0.5*rho*V^2

--- FIXED AIRCRAFT PARAMETERS ---
  m    = 1200  kg       S    = 18   m^2
  cbar = 1.5   m        Iy   = 2100 kg*m^2
  CLa  = 4.8             CLde = 0.36
  Cma  = -0.5            Cmq  = -8.0   Cmde = -1.1

--- FLIGHT ENVELOPE (keep symbolic, sweep numerically) ---
  V   :  40 – 120 m/s    (airspeed range)
  rho :  1.225 – 0.66 kg/m^3  (sea level to 5 km altitude)

Requirements
Use Symbolic Math Toolbox to:
1. Define A and B with Zalpha, Malpha, Mq, Zde, Mde as symbolic; compute the symbolic characteristic polynomial det(sI - A) and display it
2. Derive symbolic expressions for the short‑period natural frequency wn and damping ratio zeta in terms of the dimensional derivatives; display both
3. Substitute the parameter‑to‑derivative relationships so wn and zeta are symbolic functions of V and rho only (all other aircraft parameters numeric); simplify and display
4. Use matlabFunction() to export wn(V,rho) and zeta(V,rho) as numeric functions; evaluate on a grid of V and rho values spanning the flight envelope; generate surface or contour plots of wn and zeta vs. airspeed and altitude
5. Overlay MIL‑SPEC boundaries (Level 1 handling qualities: 0.35 ≤ zeta ≤ 1.3, wn ≥ 0.6 rad/s) on the damping‑ratio contour plot to show where the open‑loop aircraft meets or violates requirements
6. Save the symbolic wn, zeta expressions, the numeric matlabFunction handles, and the evaluation grid data to a MAT‑file
```

### Example 2: Black-Scholes Option Pricing

```
Create a concise plain‑text MATLAB Live Code file (.m) showing how to derive option Greeks symbolically from the Black‑Scholes formula, then export numeric functions for real‑time portfolio risk evaluation.

--- MODEL DEFINITION ---
Black‑Scholes European call price:
  C(S,K,r,sigma,T) = S*N(d1) - K*exp(-r*T)*N(d2)
where:
  d1 = (ln(S/K) + (r + sigma^2/2)*T) / (sigma*sqrt(T))
  d2 = d1 - sigma*sqrt(T)
  N(x) = standard normal CDF = (1/2)*(1 + erf(x/sqrt(2)))
Variables (keep symbolic throughout derivation):
  S      spot price
  K      strike price
  r      risk‑free rate
  sigma  volatility
  T      time to expiry
Greeks to derive:
  Delta  = dC/dS
  Gamma  = d2C/dS2
  Theta  = -dC/dT      (note sign convention: dC/dT then negate)
  Vega   = dC/d(sigma)
  Rho    = dC/dr

--- REFERENCE PARAMETER VALUES (for numeric verification) ---
  S     = 100          (spot price)
  K     = 105          (strike, slightly OTM)
  r     = 0.05         (5% risk‑free rate)
  sigma = 0.20         (20% annual volatility)
  T     = 0.5   years  (6 months to expiry)

Requirements
Use Symbolic Math Toolbox to:
1. Define C(S,K,r,sigma,T) symbolically using sym and erf for the normal CDF; display the full symbolic call price expression
2. Compute all five Greeks using diff(); simplify each; display the closed‑form symbolic expressions
3. Verify by substituting the reference parameter values and comparing against finite‑difference approximations (bump S by ±0.01 for Delta, etc.); display both values side‑by‑side
4. Use matlabFunction() with the 'File' option to export C and all five Greeks as a single optimized .m function file (e.g., bsGreeks.m) that accepts (S,K,r,sigma,T) and returns [C, Delta, Gamma, Theta, Vega, Rho]; verify the exported function matches the symbolic subs() results at the reference point
5. Sweep S from 80 to 120 and sigma from 0.10 to 0.40; generate surface plots of Delta(S,sigma) and Gamma(S,sigma); generate a line plot of Theta vs. T (0.01 to 1 year) at the reference S, K, r, sigma to show time‑decay acceleration near expiry
6. Save the symbolic Greek expressions and the numeric function handles to a MAT‑file
```

## Common Patterns

```matlab
% Pattern 1: Characteristic polynomial
charpoly = det(A - sym('lambda')*eye(2))

% Pattern 2: Generating MATLAB Functions from Symbolic Expressions

syms x y z

% Function handle (anonymous function)
f = x^2 + y^2 + z^2;
fh = matlabFunction(f);
% Creates: @(x,y,z) x.^2 + y.^2 + z.^2

% Control variable order
fh = matlabFunction(f, 'Vars', {x, y, z});

% Group variables into vectors
fh = matlabFunction(f, 'Vars', {[x, y, z]});
% Creates: @(in1) in1(:,1).^2 + in1(:,2).^2 + in1(:,3).^2

% Write to file instead of returning handle
% NOTE: Optimization is automatic when writing to file.
% Do NOT use 'Optimize', true without 'File' — it will error.
% Be sure to either include all free variables in `'Vars'`, substitute constants BEFORE calling `matlabFunction`, or group parameters into a vector input
matlabFunction(f, 'File', 'myFunction', 'Vars', {x, y, z});
% Creates myFunction.m on disk with optimized intermediate variables

% Multiple outputs
syms a b
expr1 = a + b;
expr2 = a * b;
fh = matlabFunction(expr1, expr2, 'Vars', {a, b});
% Creates: @(a,b) deal(a+b, a.*b)
```

## Related Prompts

- [Linearize Nonlinear Dynamics](./linearize.md)
- [Symbolic Transforms](./transforms.md)
- [Diagram to Symbolic Model](./diagram-to-model.md)
- [Generate Plain Text Live Scripts](../live-scripts-documentation/generate-plain-text-live-script.md)

## References

- [MATLAB Documentation: Simplify Symbolic Expressions](https://www.mathworks.com/help/symbolic/simplify-symbolic-expressions.html)
- [MATLAB Documentation: Generate MATLAB Functions from Symbolic Expressions](https://www.mathworks.com/help/symbolic/generate-matlab-functions.html)
- [MATLAB Documentation: The Black–Scholes Formula for Call Option Price](https://www.mathworks.com/help/symbolic/the-black-scholes-formula-for-call-option-price.html)
- [MATLAB Documentation: Live Code File Format (.m)](https://www.mathworks.com/help/matlab/matlab_prog/plain-text-file-format-for-live-scripts.html)
---