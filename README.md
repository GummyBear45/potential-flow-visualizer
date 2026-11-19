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

## Method

Everything is derived from the complex velocity

$$u - iv = \frac{dW}{dz} = U + \sum_k \frac{A_k}{z - z_k} + \sum_k \frac{B_k}{(z - z_k)^2}$$

with $A_k = (m_k - i\Gamma_k)/2\pi$ for sources/sinks and vortices and $B_k = -\kappa_k e^{i\alpha_k}/2\pi$ for doublets. The velocity field is single-valued, unlike the stream function and the potential, so streamlines are integrated directly from it. Stagnation points and dividing streamlines are computed from the same expression.

The results were validated against the analytical solutions for the Rankine half-body, the Rankine oval and the flow around a cylinder.
