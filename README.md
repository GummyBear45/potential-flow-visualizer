# Interactive Potential Flow Visualizer (MATLAB)

A real-time 2D potential-flow sandbox. Sources, sinks, point vortices and doublets are superposed on a uniform stream and can be added, dragged, strengthened and removed with the mouse. Velocity vectors, streamlines, stagnation points and animated tracer particles update live.

## Run

```matlab
przeplywer
```

Requires MATLAB (no toolboxes).

## Controls

| Input | Action |
|---|---|
| Left-drag | Move an element |
| Scroll | Change the strength of the element under the cursor |
| Shift + scroll | Rotate a doublet axis |
| Right-click | Open the properties panel: position, strength, doublet angle, blinking, delete |
| `Delete` | Remove the selected element |
| `s` / `k` | Add a source / sink at the cursor |
| `v` / `c` | Add a counter-clockwise / clockwise vortex |
| `d` | Add a doublet |

## Blinking elements and chaotic advection

Any element can be set to jump between two positions every half of the period *Blink T* (checkbox in the right-click panel; the next position is drawn as a hollow marker). Between jumps the flow is steady, but the time-periodic switching makes particle paths chaotic.

The **Blinking vortex demo** button sets up Aref's blinking vortex (J. Fluid Mech. 143, 1984): no freestream, one vortex with $\Gamma = 3$ alternating between $(-0.6, 0)$ and $(0.6, 0)$, and a disc of two-colour dye that is never recycled, so the stretching and folding of the dye builds up on screen. Clear returns to normal speed-coloured particles.

## Airfoil app

```matlab
profil
```

A second app for airfoils in potential flow, with two modes and a shared pressure plot.

**Joukowski airfoil (exact).** The flow around a circle in an auxiliary plane is the same freestream + doublet + vortex that `przeplywer` builds; it is shown in the bottom-right panel. The map $z = \zeta + b^2/\zeta$ turns the circle into an airfoil. Sliders move the circle centre (left for thickness, up for camber) and set the angle of attack. The Kutta condition fixes the circulation, so the flow leaves the cusped trailing edge smoothly. Streamlines and particles are computed in the circle plane and mapped.

**Thin airfoil (source/vortex distributions).** Sources and vortices are spread along the chord as constant-strength panels. Drag the control points of the source density (red) and vortex density (green): the thickness follows from the sources and the camber line from the vortices' downwash, and the airfoil redraws live. The last vortex point is locked at zero (Kutta condition). Presets give a symmetric NACA-0012-like section, a 4 % parabolic camber at its ideal angle, and a plate at 5°. **Keep body closed** (on by default) shifts all the other source points, except the leading-edge one, while you drag, so the net source strength stays zero and the trailing edge stays closed. With it off, **Close body** closes the body once. The info box warns when the thickness turns negative. The plate preset shows about 4° between leading and trailing edge, not 5°: finite control points cannot hold the leading-edge singularity of the plate's vortex sheet, and cutting it off droops the nose.

**Pressure coefficient.** The right-hand plot shows $C_p$ along the chord, negative up. In thin airfoil theory $C_p$ is linear in the disturbance, so it splits exactly into a thickness part (the same on both surfaces), a camber part and an angle-of-attack part (opposite signs on the two surfaces). In Joukowski mode the exact $C_p$ is drawn over the thin-theory curves for the same shape, and the lift coefficients of both are listed, so the effect of thickness on lift is visible directly.

## Code layout

| File | Contents |
|---|---|
| `przeplywer.m` | free-element sandbox |
| `profil.m` | airfoil app |
| `+pf/` | shared core: pole velocity, streamline tracer, capped arrows, particle step, stagnation points, figure setup |

## Method

Everything is derived from the complex velocity

$$u - iv = \frac{dW}{dz} = U + \sum_k \frac{A_k}{z - z_k} + \sum_k \frac{B_k}{(z - z_k)^2}$$

with $A_k = (m_k - i\Gamma_k)/2\pi$ for sources/sinks and vortices and $B_k = -\kappa_k e^{i\alpha_k}/2\pi$ for doublets. The velocity field is single-valued, unlike the stream function and the potential, so streamlines are integrated directly from it. Stagnation points and dividing streamlines are computed from the same expression.

The results were validated against the analytical solutions for the Rankine half-body, the Rankine oval and the flow around a cylinder.

The airfoil app was checked against known results: a flat plate gives $C_L = 2\pi \sin\alpha$ exactly, the contour integral of the velocity around a Joukowski airfoil equals the Kutta circulation, the velocity stays finite at the trailing edge, and a symmetric section at zero angle of attack gives zero lift. The thin-theory presets give $C_L = 0$ (symmetric) and $C_L = 0.50$ (4 % parabolic camber, theory $0.503$).
