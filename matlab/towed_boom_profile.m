function profile = towed_boom_profile(alpha_deg, T0, tau, nphi)
% TOGED_BOOM_PROFILE Compute the shape of a towed boom of logs using
% elliptic-integral formulas from Newman (as summarized in the ticket).
%
% profile = towed_boom_profile(alpha_deg, T0, tau, nphi)
%   alpha_deg : half of the included angle at the tow point, in degrees
%               0 < alpha_deg < ~41
%   T0        : constant boom tension (units of force)
%   tau       : skin-friction (units of force per unit area)
%   nphi      : number of samples along the profile (default 400)
%
% Returns a struct containing fields:
%   alpha, m, A, L, K, Ecomplete
%   phi, x, z, s               (dimensional)
%   x_scaled, z_scaled, s_scaled (scaled by sqrt(A))
%
% Key formulas implemented:
%   m = (1 + sin(alpha))/2
%   A = 2*T0*cos(alpha)/tau = 4*(T0/tau)*sqrt(m*(1-m))
%   s(phi) = sqrt(T0/tau)*(K(m) - F(phi|m))
%   x(phi) = sqrt(T0/tau)*(2E(phi|m) - F(phi|m))
%   z(phi)/sqrt(A) = sqrt(m)*cos(phi) / (m*(1-m))^(1/4)
%   L/sqrt(A) = 2/((2m-1)*(m*(1-m))^(1/4)) * (E(π/2|m) - (1-m)K(m))
%
% The full boom consists of two symmetric branches about z=0. The arrays
% here give one branch (z >= 0); mirror for the other branch if desired.

    if nargin < 4 || isempty(nphi)
        nphi = 400;
    end
    % Basic input checks (no 'arguments' block for wider compatibility)
    if ~isscalar(alpha_deg) || ~isfinite(alpha_deg)
        error('alpha_deg must be a finite scalar.');
    end
    if ~isscalar(T0) || ~isfinite(T0) || T0 <= 0
        error('T0 must be a positive finite scalar.');
    end
    if ~isscalar(tau) || ~isfinite(tau) || tau <= 0
        error('tau must be a positive finite scalar.');
    end
    if ~isscalar(nphi) || ~isfinite(nphi) || nphi <= 0 || nphi ~= round(nphi)
        error('nphi must be a positive integer scalar.');
    end

    alpha = alpha_deg*pi/180;
    if alpha <= 0 || alpha >= 41.0*pi/180
        warning('alpha_deg should satisfy 0 < alpha_deg < ~41. Using provided value: %.3f deg', alpha_deg);
    end

    % Parameter m in (1/2, 1)
    m = (1 + sin(alpha))/2;

    % Sample parameter phi in [0, pi/2]
    phi = linspace(0, pi/2, nphi);

    % Complete elliptic integrals K(m), E(m)
    [Kc, Ec] = ellipke_strict(m);

    % Incomplete elliptic integrals F(phi|m), E(phi|m)
    Fphi = ellipticF_incomplete(phi, m);
    Ephi = ellipticE_incomplete(phi, m);

    % Scale factor
    rt = sqrt(T0/tau);

    % Formulas (1.56), (1.58)
    s = rt * (Kc - Fphi);
    x = rt * (2.*Ephi - Fphi);

    % Area from overall equilibrium and from (1.59)
    A_from_equil = 2*T0*cos(alpha)/tau;
    % Also consistent with (1.59): A = 4*(T0/tau)*sqrt(m*(1-m))
    A_from_159 = 4*(T0/tau)*sqrt(m*(1-m));

    % Prefer the equilibrium expression for numerical stability
    A = A_from_equil;

    % z from (1.60); using an equivalent simplified expression to
    % improve numerical conditioning near m -> 1:
    % z = 2*sqrt(T0/tau)*sqrt(m)*cos(phi)
    z = 2*rt*sqrt(m)*cos(phi);

    % Scaled coordinates by sqrt(A)
    rootA = sqrt(A);
    x_scaled = x / rootA;
    s_scaled = s / rootA;
    z_scaled = z / rootA;

    % Verify consistency of z_scaled with (1.60)
    % z_scaled_expected = sqrt(m).*cos(phi) ./ ( (m*(1-m))^(1/4) ); %#ok<NASGU>

    % Total boom length from (1.61)
    denom = (2*m - 1) * (m*(1-m))^(1/4);
    if abs(denom) < 1e-12
        L_scaled = NaN; % too close to singular lower bound (alpha->0)
    else
        L_scaled = 2 * (Ec - (1 - m)*Kc) / denom;
    end
    L = rootA * L_scaled;

    % Populate output struct
    profile = struct();
    profile.alpha = alpha;
    profile.alpha_deg = alpha_deg;
    profile.m = m;
    profile.A = A;
    profile.A_from_equil = A_from_equil;
    profile.A_from_159 = A_from_159;
    profile.L = L;
    profile.L_scaled = L_scaled;
    profile.K = Kc;
    profile.Ecomplete = Ec;

    profile.phi = phi;
    profile.x = x;
    profile.z = z;
    profile.s = s;

    profile.x_scaled = x_scaled;
    profile.z_scaled = z_scaled;
    profile.s_scaled = s_scaled;

end

function [Kc, Ec] = ellipke_strict(m)
% Wrapper around ellipke that validates m in [0,1)
    if m < 0 || m >= 1
        error('Parameter m must satisfy 0 <= m < 1');
    end
    try
        [Kc, Ec] = ellipke(m);
    catch
        % Octave compatibility
        [Kc, Ec] = ellipke(m);
    end
end

function F = ellipticF_incomplete(phi, m)
% Compute incomplete elliptic integral of the first kind F(phi|m)
    if exist('ellipticF', 'file') == 2 || exist('ellipticF', 'builtin')
        F = ellipticF(phi, m);
    else
        % Fallback: numeric quadrature
        F = arrayfun(@(p) integral(@(t) 1./sqrt(1 - m.*sin(t).^2), 0, p, 'RelTol',1e-9,'AbsTol',1e-12), phi);
    end
end

function E = ellipticE_incomplete(phi, m)
% Compute incomplete elliptic integral of the second kind E(phi|m)
    if exist('ellipticE', 'file') == 2 || exist('ellipticE', 'builtin')
        E = ellipticE(phi, m);
    else
        % Fallback: numeric quadrature
        E = arrayfun(@(p) integral(@(t) sqrt(1 - m.*sin(t).^2), 0, p, 'RelTol',1e-9,'AbsTol',1e-12), phi);
    end
end
