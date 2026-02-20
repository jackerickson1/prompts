---
title: Convert Diagrams to Symbolic Models
description: Extract components and topology from an uploaded image, build a symbolic model, and export numeric functions or transfer functions for downstream analysis.
tags: [matlab, symbolic-math, image-input, matlabFunction, live-script]
release: All releases / Requires Symbolic Math Toolbox
notes:
---

# Convert Diagrams to Symbolic Models

Uses an AI assistant's vision capability to interpret an uploaded diagram (circuit schematic, robot arm sketch, stress-strain curve, etc.), then translates the identified structure into symbolic equations with Symbolic Math Toolbox. Useful for quickly turning hand-drawn sketches, textbook figures, or datasheet plots into executable MATLAB models without manual equation entry.

## Metadata

- **Tags:** `matlab` `symbolic-math` `image-input` `matlabFunction` `live-script`
- **MATLAB Release:** All releases
- **Required Toolboxes:** Symbolic Math Toolbox

## The Prompt

```text
I have uploaded an image of [A CIRCUIT SCHEMATIC | A ROBOT ARM DIAGRAM | A STRESS‑STRAIN CURVE | A BLOCK DIAGRAM | OTHER].
Create a concise plain‑text MATLAB Live Code file (.m) that builds a symbolic model from the image and produces [ENGINEERING DELIVERABLE: numeric tf, callable function file, fitted model, etc.].

--- FROM THE IMAGE ---
1. Identify all components, variables, or features visible in the image. List labels, values, and units where readable.
2. Identify the topology, structure, or behavior type (e.g., filter topology, kinematic chain, material response class).
3. For any values not visible in the image, declare symbolic parameters so the result stays general.
List all identifications in a comment block at the top of the script.

--- EXAMPLE (use if no image is provided) ---
  [PROVIDE A CONCRETE FALLBACK: component list, parameter table, or sample data points so the prompt is self‑contained and testable without an image]

Requirements
Use Symbolic Math Toolbox to:
1. Translate the identified structure into symbolic equations ([KVL/KCL | transformation matrices | constitutive model | transfer function block algebra]); display the symbolic model
2. Derive key engineering quantities symbolically ([corner frequencies, Q | Jacobian, singularities | tangent modulus, yield criterion]); display
3. Substitute any numeric values read from the image; where a downstream numeric object is appropriate ([tf | ss | fitted curve]), build it and display
4. Generate verification or summary plots ([Bode | workspace reachability | model overlay on data])
5. Use matlabFunction() with the 'File' option to export the symbolic model as a callable .m function file for downstream use; verify the exported function at a test point
6. Save symbolic expressions, numeric objects, and identified parameters to a MAT‑file
```

## Usage Tips

1. **Providing a diagram:** Try to ensure that the your image is high resolution, and add to the prompt any variables and values not captured clearly in the image. And be mindful of uploading confidential diagrams to public AI models.
2. **Live Script Formatting:** Symbolic equations display best as Live Script output. Starting with R2025a, MATLAB supports plain text Live Code formatting. For best results, use the [MATLAB Plain Text Live Script Generator](https://github.com/matlab/skills/blob/main/skills/matlab-live-script/SKILL.md) skill. Alternatively, you can include these instructions:
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
3. **Symbolic Math Toolbox skill:** For best results, load the [Symbolic Math Toolbox skill](https://github.com/matlab/skills/blob/main/skills/symbolic-math-toolbox/SKILL.md)
4. **Display equations and generate plots:** For each step, generate outputs so you can verify the generated code.

## Example Usage

### Example 1: RLC Circuit
```
I have uploaded an image of a circuit schematic.
Create a concise plain‑text MATLAB Live Code file (.m) that derives the symbolic transfer function and produces a numeric model for frequency analysis.

--- FROM THE IMAGE ---
1. Identify every component (R, L, C, op‑amp, source, etc.) and its label/value. List them in a comment block at the top of the script.
2. Identify the input node(s) and output node(s).
3. Identify the circuit topology (voltage divider, feedback network, filter stage, etc.).
If any component values are missing from the image, declare them as symbolic parameters so the result stays general.

--- NUMERICAL PARAMETERS ---
Extract from the image

Requirements
Use Symbolic Math Toolbox to:
1. Define symbolic variables for voltages, currents, inductances, capacitances, and the Laplace variable s
2. Write KVL/KCL equations in the Laplace domain to create the system of equations. Display the system of equations.
3. Solve the system of equations for the variables I1, I2, V1, and V2 and display the result.
4. Derive the transfer function H(s) = Vout(s)/Vin(s)
5. Substitute the numerical parameters into the transfer function and use variable precision to compute the transfer function to 5 digits of accuracy.
6. Extract the numerator and denominator of the transfer function, convert them to numeric vectors of polynomial coefficients, and use them to construct a tf object. Cancel any pole-zero pairs, and display the tf.
7. Plot the pole-zero map. Create a Bode plot for the magnitude and phase response for frequencies linspace(1e11,1e15,4000).
```

### Example 2: Robotic Arm Forward Kinematics

```
I have uploaded a diagram of a serial robot arm.

Create a concise plain‑text MATLAB Live Code file (.m) that derives symbolic forward kinematics and the geometric Jacobian, then exports fast numeric functions for use in a control loop.

--- FROM THE IMAGE ---
1. Identify the number of joints, joint types (revolute or prismatic), and label each joint variable (q1, q2, …).
2. Identify link lengths, offsets, and any fixed angles visible in the diagram. Assign symbolic parameters for any dimensions not labeled.
3. Assign DH parameters for each link, or use the geometry directly if the arm is planar.
List all identifications in a comment block at the top of the script.

Requirements

Use Symbolic Math Toolbox to:
1. Build symbolic homogeneous transformation matrices for each link using DH convention. Multiply to get the overall forward kinematics T_0n, then extract the end‑effector position and orientation. Display T_0n.
2. Compute the geometric Jacobian J(q) = dp/dq, then simplify and display it.
3. Evaluate the Jacobian determinant symbolically to identify singular configurations. Solve det(J) = 0 if feasible and display the result.
4. Substitute numeric link dimensions from the image. Evaluate forward kinematics and Jacobian at all joints equal to 0. Display numeric end‑effector pose and Jacobian.
5. Use matlabFunction() to export the forward kinematics and Jacobian as optimized .m functions that accept a joint vector and return numeric results. Verify by calling them at the same test configuration.
```

### Example 3: Stress-Strain Curve to Symbolic Constitutive Model

```
I have uploaded an image of a stress‑strain curve.
Create a concise plain‑text MATLAB Live Code file (.m) that proposes a symbolic constitutive model, derives engineering quantities, and exports
a numeric function for curve fitting.

--- FROM THE IMAGE ---
1. Identify the type of material behavior visible in the curve:
   linear elastic, elastoplastic with hardening, hyperelastic (rubber), viscoelastic, or other. Note any labeled axes, units, and data points.
2. Identify approximate key values from the plot - elastic modulus, yield stress, ultimate stress, and failure strain
3. Propose an appropriate symbolic constitutive model, for instance:
   - Linear + power‑law hardening
   - Ramberg‑Osgood
   - Mooney‑Rivlin (hyperelastic)
   - Or another model that best matches the observed curve shape
List the identification and chosen model in a comment block at the top of the script.

Requirements
Use Symbolic Math Toolbox to:
1. Define the chosen constitutive model symbolically with material parameters as syms, then display the stress‑strain relationship sigma(eps).
2. Derive the symbolic tangent modulus, then simplify and display it.
3. Derive and display symbolic equations for: 
   specific strain energy
   onset of necking condition
4. Use matlabFunction() to convert the symbolic sigma(eps) and E_t(eps) into numeric functions that accept a strain vector and parameter values.
5. Generate representative plots: overlay 2–3 curves with different parameter values to show model sensitivity (e.g., vary the hardening exponent n); plot the tangent modulus vs. strain alongside.
```

## Common Patterns

```matlab
% Pattern 1: Direct Symbolic Transfer Function

syms s K wn zeta

% Standard second-order system
G = K * wn^2 / (s^2 + 2*zeta*wn*s + wn^2);

% Substitute numeric values
G_numeric = subs(G, [K, wn, zeta], [1, 10, 0.5]);

% Convert to tf object
[n, d] = numden(G_numeric);
sys = tf(sym2poly(n), sym2poly(d));

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

- [Symbolic Transforms](./transforms.md)
- [Transforms and Spectral Methods](./transforms.md)
- [Symbolic Parameters in Numerical Workflows](./parametric.md)
- [Generate Plain Text Live Scripts](../live-scripts-documentation/generate-plain-text-live-script.md)

## References

- [MATLAB Documentation: Simplify Symbolic Expressions](https://www.mathworks.com/help/symbolic/simplify-symbolic-expressions.html)
- [MATLAB Documentation: Choose an Approach for Solving Equations Using solve Function](https://www.mathworks.com/help/symbolic/troubleshoot-equation-solutions-from-solve-function.html)
- [MATLAB Documentation: Analyze Transfer Function of T-Coil Circuit](https://www.mathworks.com/help/symbolic/analyze-transfer-function-of-t-coil.html)
- [MATLAB Documentation: Derive and Apply Inverse Kinematics to Two-Link Robot Arm](https://www.mathworks.com/help/symbolic/derive-and-apply-inverse-kinematics-to-robot-arm.html)
- [MATLAB Documentation: Live Code File Format (.m)](https://www.mathworks.com/help/matlab/matlab_prog/plain-text-file-format-for-live-scripts.html)
---