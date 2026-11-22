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

**Thin airfoil (source/vortex distributions).** Sources and vortices are spread along the chord as constant-strength panels. Drag the control points of the source density (red) and vortex density (green): the thickness follows from the sources and the camber line from the vortices' downwash, and the airfoil redraws live. The last vortex point is locked at zero (Kutta condition). Because thin airfoil theory is linear, the three sliders act directly: **thickness scale** multiplies the source density (and so the thickness), **camber scale** multiplies the vortex density (camber and its lift), and **angle of attack** adds the exact flat-plate vortex sheet $\gamma = 2U\alpha\sqrt{(1-x)/x}$, whose downwash is a uniform $-U\alpha$, so the airfoil turns nose-up in the horizontal freestream. A fourth slider (bottom, shared with the trailing-edge angle of the Joukowski mode) sets the **nose radius** $r$: it adds the exact round-nose thickness $y_N = \sqrt{2rx}\,(1-x)$, which closes at the trailing edge by itself, and its surface speed is computed in closed form because constant-strength panels cannot hold its $1/\sqrt{x}$ source density. Presets give an exact NACA 0012 (round-nose term $r/c = 0.0159$ plus a smooth remainder in the editor), a 4 % parabolic camber on the same thickness at its ideal angle, and a sharp flat plate at 5°. **Keep body closed** (on by default) shifts all the other source points, except the leading-edge one, while you drag, so the net source strength stays zero and the trailing edge stays closed. The info box warns when the thickness turns negative.

**Trailing-edge angle.** The Joukowski airfoil has a cusped trailing edge. The *TE angle* slider replaces the map by the Karman-Trefftz map $(z - nb)/(z + nb) = ((\zeta - b)/(\zeta + b))^n$ with $n = 2 - \tau/\pi$, which gives a trailing edge of angle $\tau$ and keeps the solution exact ($\tau = 0$ is the Joukowski airfoil). This plays the role of XFOIL's TGAP; an open (blunt) trailing edge itself cannot be made with an exact conformal map, because the potential flow around a blunt edge is singular.

**Boundary layer (optional, both modes).** Tick *Boundary layer* and enter a Reynolds number. The surface speed drives Thwaites' laminar method from the stagnation point, Michel's criterion (in its full form, which cannot fire next to the stagnation point) or laminar separation switches to Head's turbulent entrainment method, and Squire-Young gives the profile drag $c_d$ and $L/D$. The displacement thickness is drawn five times enlarged, triangles mark transition and stars separation. It is stopped at $x/c = 0.99$ (an inviscid finite-angle edge is a stagnation point, a real wake is not). In Joukowski mode the surface speed is exact, the boundary layer does not feed back on the outer flow, and drag is not reported once the flow separates; thin-airfoil mode couples the two (next paragraph). In thin-airfoil mode linear theory has no stagnation point and an infinite suction peak at the nose, so the speed uses the uniformly valid leading-edge correction of Lighthill and Van Dyke,

$$\frac{q}{U} = \left(1 + \frac{u_t}{U} \pm \frac{\gamma_c}{2U}\right)\sqrt{\frac{x}{x + r/2}} \pm \alpha\sqrt{\frac{1 - x}{x + r/2}},$$

with $u_t$ the thickness-induced speed, $\gamma_c$ the editor (camber) vortex density, $\alpha$ the flat-plate part and $r$ the nose radius; the stagnation point is where this speed changes sign. A sharp nose ($r = 0$) at incidence therefore separates at once, as it should.

**Viscous-inviscid interaction (thin-airfoil mode).** Past separation the pressure is no longer the potential-flow pressure, so the boundary layer is coupled back into the outer flow:

- the displacement thickness $\delta^\ast$ of both surfaces and of a wake 1.5 chords long acts on the outer flow as extra sources $U\,\mathrm{d}\delta^\ast/\mathrm{d}x$ on the chord line and the wake axis, and as an extra camber $(\delta^\ast_u - \delta^\ast_l)/2$, re-solved with Glauert's series, which lowers the lift;
- the turbulent layer is marched directly (outer speed given) while attached and, once its shape factor passes 2, in inverse mode ($\delta^\ast$ given, speed solved), which stays regular through separation where the direct method is singular;
- $\delta^\ast$ is corrected with Carter's semi-inverse rule $\delta^\ast \leftarrow \delta^\ast\,[1 + \omega(U_{e,BL}/U_{e,outer} - 1)]$, with a local relaxation $\omega \le \Delta s/(\pi\delta^\ast)$ (the outer response to a wiggle grows with its wavenumber), until both sides agree to 0.5 %;
- the wake is marched with Head's method without wall friction, and the drag is Squire-Young at its far end, so $c_d$ stays defined with separation on the airfoil;
- the solution is warm-started from the last one; if a large jump fails, the angle of attack is stepped up from the last good angle (continuation). The info box says whether the coupling converged.

The last 4 % of chord is left out of the coupling: at a finite-angle trailing edge linear theory slows the flow to zero, and the iteration would try to fill that with displacement. The viscous surface pressure (leading-edge corrected, with displacement) replaces the "exact" curves in the $C_p$ plot. NACA 0012 at $Re = 10^6$: $c_d$ = 0.0073, 0.0080, 0.0100, 0.0135 and 0.0191 at 0°, 4°, 8°, 12° and 16°, with viscous $C_L$ 5–8 % below the inviscid value and separation appearing on the upper surface from about 16° (x/c 0.95, moving to 0.83 at 18°). Measured and XFOIL values are close up to about 12°; beyond that Head's integral method under-predicts the separated region, so stall and post-stall drag come out too mild.

Checks: a nearly flat Joukowski section at $Re = 2 \cdot 10^5$ gives $c_d = 0.00616$ against the laminar flat-plate (Blasius) value $0.00594$. For the NACA 0012 in thin-airfoil mode the stagnation point moves to $x/c$ = 0.0008, 0.0048 and 0.0124 on the lower surface at 2°, 5° and 8°, the suction peaks are $C_p$ = −0.63, −1.50 and −3.04, transition lies at $x/c \approx 0.5$ at 0° and at 0.17 (upper) / 0.83 (lower) at 5° for $Re = 10^6$, and $c_d$ at 0° is 0.0063 at $Re = 3 \cdot 10^6$ and 0.0058 at $6 \cdot 10^6$, against about 0.0058 and 0.0054 measured for the smooth section.

**Pressure coefficient.** The right-hand plot shows $C_p$ along the chord, negative up. In thin airfoil theory $C_p$ is linear in the disturbance, so it splits exactly into a thickness part (the same on both surfaces), a camber part and an angle-of-attack part (opposite signs on the two surfaces). In Joukowski mode the exact $C_p$ is drawn over the thin-theory curves for the same shape, and the lift coefficients of both are listed, so the effect of thickness on lift is visible directly.

## Code layout

| File | Contents |
|---|---|
| `przeplywer.m` | free-element sandbox |
| `profil.m` | airfoil app |
| `+pf/` | shared core: pole velocity, streamline tracer, capped arrows, particle step, stagnation points, figure setup, timer cleanup |

## Method

Everything is derived from the complex velocity

$$u - iv = \frac{dW}{dz} = U + \sum_k \frac{A_k}{z - z_k} + \sum_k \frac{B_k}{(z - z_k)^2}$$

with $A_k = (m_k - i\Gamma_k)/2\pi$ for sources/sinks and vortices and $B_k = -\kappa_k e^{i\alpha_k}/2\pi$ for doublets. The velocity field is single-valued, unlike the stream function and the potential, so streamlines are integrated directly from it. Stagnation points and dividing streamlines are computed from the same expression.

The results were validated against the analytical solutions for the Rankine half-body, the Rankine oval and the flow around a cylinder.

The airfoil app was checked against known results: a flat plate gives $C_L = 2\pi \sin\alpha$ exactly, the contour integral of the velocity around a Joukowski airfoil equals the Kutta circulation, the velocity stays finite at the trailing edge, and a symmetric section at zero angle of attack gives zero lift. The thin-theory presets give $C_L = 0$ (symmetric), $C_L = 0.50$ (4 % parabolic camber, theory $0.503$) and an angle of 4.99° for the 5° plate; doubling either scale slider doubles the thickness or the camber and its lift exactly.
