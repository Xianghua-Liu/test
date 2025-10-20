# MATLAB reproduction: Shape of a Towed Boom of Logs (Newman)

This repository contains a minimal MATLAB implementation that reproduces the elliptic-integral solution for the shape of a towed boom of logs, following the derivation summarized in the ticket. The implementation evaluates the profile coordinates, arc length, and total boom length using incomplete and complete elliptic integrals.

Key equations used (numbering from the summary):
- m = (1 + sin α)/2, with 1/2 < m < 1
- s(φ) = (T₀/τ)^{1/2} { K(m) − F(φ|m) } (eq. 1.56)
- x(φ) = (T₀/τ)^{1/2} { 2E(φ|m) − F(φ|m) } (eq. 1.58)
- A = 2 T₀ cos α / τ = 4 (T₀/τ) [m(1−m)]^{1/2} (eqs. Overall, 1.59)
- z/A^{1/2} = m^{1/2} cos φ / [m(1−m)]^{1/4} (eq. 1.60)
- L/A^{1/2} = 2/[(2m−1)[m(1−m)]^{1/4}] { E(π/2|m) − (1−m)K(m) } (eq. 1.61)

The code also captures the symmetry of the boom about the z=0 axis by plotting both the upper and mirrored lower branches.

## Contents

- `matlab/towed_boom_profile.m` — function computing the profile given α, T₀, τ
- `matlab/demo_towed_boom_profiles.m` — script that plots several dimensionless profiles for different included angles

## Requirements

- MATLAB (recommended) with elliptic integral functions `ellipticF`, `ellipticE`, and `ellipke`
  - If `ellipticF`/`ellipticE` are unavailable, the code falls back to accurate numeric quadrature for the incomplete integrals.
- GNU Octave should also work, leveraging the numeric quadrature fallback for incomplete integrals and `ellipke` for complete integrals.

## Quick start

1. Open MATLAB (or Octave) and add the `matlab` folder to your path:

   addpath(genpath('matlab'))

2. Run the demo to plot several profiles (dimensionless scaling x/√A vs z/√A):

   demo_towed_boom_profiles

You should see a figure with several profiles corresponding to different values of the included angle (2α). The script also prints some representative scalar outputs for one case.

## Programmatic use

Call the core function directly to compute a profile structure with fields for coordinates, scaled coordinates, area A, total length L, and more:

   % Inputs: alpha_deg (half-included angle in degrees), T0, tau
   prof = towed_boom_profile(20, 1.0, 1.0, 600);

   % Dimensional coordinates along the curve parameterized by φ ∈ [0, π/2]
   x = prof.x;     % along-tow coordinate
   z = prof.z;     % transverse coordinate (upper branch)
   s = prof.s;     % arc length from φ to the apex

   % Dimensionless, scaled by √A
   xs = prof.x_scaled;
   zs = prof.z_scaled;
   ss = prof.s_scaled;

   % Scalars
   A  = prof.A;
   L  = prof.L;
   m  = prof.m;

Note: The shape family is governed by α (or equivalently by m). The ratio (T₀/τ) sets the overall length/area scale via A ∝ (T₀/τ).

## Notes

- Valid range for α is 0 < α ≲ 41° (so that 1/2 < m < 1). As α approaches its upper bound, m → 1 and the profile becomes very full; as α → 0, the expression for total length (eq. 1.61) becomes singular.
- The demo uses T₀ = τ = 1. Changing (T₀/τ) rescales x, z, s, and L by (T₀/τ)^{1/2}, and rescales A linearly.

