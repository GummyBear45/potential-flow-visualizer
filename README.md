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

**Thin airfoil (source/vortex distributions).** Sources and vortices are spread along the chord as constant-strength panels; their strengths come from linear theory (on $y = 0$), and for the flow picture the panels are placed along the mean line so that the streamlines follow the airfoil drawn. Drag the control points of the source density (red) and vortex density (green): the thickness follows from the sources and the camber line from the vortices' downwash, and the airfoil redraws live. The last vortex point is locked at zero (Kutta condition). Because thin airfoil theory is linear, the three sliders act directly: **thickness scale** multiplies the source density (and so the thickness), **camber scale** multiplies the vortex density (camber and its lift), and **angle of attack** adds the exact flat-plate vortex sheet $\gamma = 2U\alpha\sqrt{(1-x)/x}$, whose downwash is a uniform $-U\alpha$, so the airfoil turns nose-up in the horizontal freestream. A fourth slider (bottom, shared with the trailing-edge angle of the Joukowski mode) sets the **nose radius** $r$: it adds the exact round-nose thickness $y_N = \sqrt{2rx}\thinspace (1-x)$, which closes at the trailing edge by itself, and its surface speed is computed in closed form because constant-strength panels cannot hold its $1/\sqrt{x}$ source density. Presets give an exact NACA 0012 (round-nose term $r/c = 0.0159$ plus a smooth remainder in the editor), a 4 % parabolic camber on the same thickness at its ideal angle, and a sharp flat plate at 5°. **Keep body closed** (on by default) shifts all the other source points, except the leading-edge one, while you drag, so the net source strength stays zero and the trailing edge stays closed. The info box warns when the thickness turns negative.

**Trailing-edge angle.** The Joukowski airfoil has a cusped trailing edge. The *TE angle* slider replaces the map by the Karman-Trefftz map $(z - nb)/(z + nb) = ((\zeta - b)/(\zeta + b))^n$ with $n = 2 - \tau/\pi$, which gives a trailing edge of angle $\tau$ and keeps the solution exact ($\tau = 0$ is the Joukowski airfoil). This plays the role of XFOIL's TGAP; an open (blunt) trailing edge itself cannot be made with an exact conformal map, because the potential flow around a blunt edge is singular.

**Boundary layer (optional, thin-airfoil mode).** Tick *Boundary layer* and enter a Reynolds number. The surface speed drives Thwaites' laminar method from the stagnation point. Transition is predicted with the $e^N$ method: the amplification $N$ of the most unstable Tollmien-Schlichting waves is integrated along the surface with the envelope correlation of Drela and Giles (the one XFOIL uses), and the layer turns turbulent at $N = 9$ (XFOIL's default $N_{crit}$, for a quiet wind tunnel). If the laminar layer separates first, it continues as a short laminar separation bubble whose shape factor climbs from 3.55, so the waves grow much faster and transition follows within a few per cent of chord, as in XFOIL; the bubble is limited to 0.1 chord. After transition Head's turbulent entrainment method takes over, and Squire-Young gives the profile drag $c_d$ and $L/D$. The displacement thickness is drawn five times enlarged, triangles mark transition and stars separation. It is stopped at $x/c = 0.99$ (an inviscid finite-angle edge is a stagnation point, a real wake is not). Joukowski mode shows the potential flow only (no boundary layer): its cusped trailing edge makes the coupling unreliable once the flow separates.  In thin-airfoil mode linear theory has no stagnation point and an infinite suction peak at the nose, so the speed uses the uniformly valid leading-edge correction of Lighthill and Van Dyke,

$$\frac{q}{U} = \left(1 + \frac{u_t}{U} \pm \frac{\gamma_c}{2U}\right)\sqrt{\frac{x}{x + r/2}} \pm \alpha\sqrt{\frac{1 - x}{x + r/2}},$$

with $u_t$ the thickness-induced speed, $\gamma_c$ the editor (camber) vortex density, $\alpha$ the flat-plate part and $r$ the nose radius; the stagnation point is where this speed changes sign. A sharp nose ($r = 0$) at incidence therefore separates at once, as it should.

For the boundary layer the two lifting terms (the $\gamma_c$ and $\alpha$ terms) are further multiplied by the local thickness factor $1 + u_t/U$, a second-order correction: the lifting flow is carried by the faster flow over the thick section. For an ellipse this is exact at any angle, because there $u_t/U$ equals the thickness ratio everywhere. It raises the NACA 0012 suction peak at 8° from $C_p = -3.04$ to $-4.39$ (XFOIL inviscid: $-4.27$) and gives the known gain in lift slope from thickness, about $1 + 0.77\thinspace t/c$. With the boundary layer on, $C_L$ is therefore the integral of the surface pressure, not the thin-theory circulation; the info box gives it with and without the displacement effect.

**Viscous-inviscid interaction (thin-airfoil mode).** Past separation the pressure is no longer the potential-flow pressure, so the boundary layer is coupled back into the outer flow:

- the displacement thickness $\delta^\ast$ of both surfaces and of a wake 1.5 chords long acts on the outer flow as extra sources $U\thinspace \mathrm{d}\delta^\ast/\mathrm{d}x$ on the chord line and the wake axis, and as an extra camber $(\delta^\ast_u - \delta^\ast_l)/2$, re-solved with Glauert's series, which lowers the lift;
- the turbulent layer is marched directly (outer speed given) while attached and, once its shape factor passes 2, in inverse mode ($\delta^\ast$ given, speed solved), which stays regular through separation where the direct method is singular;
- transition, laminar separation, the stagnation panel and the switch between the direct and inverse turbulent methods are frozen after a few iterations, so no station keeps flipping between two states;
- $\delta^\ast$ is corrected with Carter's semi-inverse rule $\delta^\ast \leftarrow \delta^\ast\thinspace [1 + \omega(U_{e,BL}/U_{e,outer} - 1)]$, with a local relaxation $\omega \le \Delta s/(\pi\delta^\ast)$ (the outer response to a wiggle grows with its wavenumber), until both sides agree to 0.5 %;
- the wake is marched with Head's method without wall friction, and the drag is Squire-Young at its far end, so $c_d$ stays defined with separation on the airfoil;
- the solution is warm-started from the last one; if a large jump fails, the angle of attack is stepped up from the last good angle (continuation). The info box says whether the coupling converged.

The last 4 % of chord is left out of the coupling: at a finite-angle trailing edge linear theory slows the flow to zero, and the iteration would try to fill that with displacement. The viscous surface pressure (leading-edge corrected, with displacement) replaces the "exact" curves in the $C_p$ plot.

**Comparison with XFOIL** (NACA 0012, $Re = 10^6$, $N_{crit} = 9$; figure `polar_vs_xfoil.png`). The coupling converges from 0° to 14°. Upper-surface transition agrees with XFOIL to within 0.01 chord up to 8° (0.69, 0.48, 0.27, 0.08 and 0.05 at 0°, 2°, 4°, 6° and 8°), and $c_d$ is within 4 % from 0° to 4° and 7 % at 5° (0.0054, 0.0059, 0.0071 and 0.0079 at 0°, 2°, 4° and 5°; XFOIL 0.0053, 0.0057, 0.0072 and 0.0085). From 6° to 12° $c_d$ is 15–28 % low (0.0100 against 0.0123 at 8°, 0.0137 against 0.0190 at 12°), and $C_L$ follows XFOIL to about 10° but does not stall: Head's integral method grows the turbulent layer too slowly in the strong pressure rise behind the suction peak and misses the trailing-edge separation that XFOIL finds from about 12°. At high angles XFOIL keeps the lower surface laminar to the trailing edge; here its bubble is closed at $x/c = 0.94$.

Checks: for the NACA 0012 in thin-airfoil mode the stagnation point (with the boundary layer on) moves to $x/c$ = 0.0011, 0.0071 and 0.0181 on the lower surface at 2°, 5° and 8°, the first-order suction peaks are $C_p$ = −0.63, −1.50 and −3.04, and $c_d$ at 0° is 0.0051 at $Re = 3 \cdot 10^6$ and 0.0052 at $6 \cdot 10^6$, against about 0.0058 and 0.0054 measured for the smooth section (with $N_{crit} = 9$ the laminar run is long, 0.51 and 0.38 chord, as for an airfoil in very quiet air).

**Polar.** *Polar...* asks for an angle range and sweeps it for the airfoil as it is set up. In thin-airfoil mode it runs the coupled boundary layer at the current Reynolds number (stepping away from 0° in both directions, so each solution starts from its neighbour) and plots $C_L$ against the angle (viscous and inviscid), the drag polar $C_L(c_d)$, $c_d$ against the angle and $L/D$; points where the coupling did not converge are drawn as open markers. In Joukowski mode it plots the exact lift against the thin-theory lift. The table is put in the MATLAB workspace as `profilPolar` and can be saved as CSV.

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
