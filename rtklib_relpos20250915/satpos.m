function [svpos, svclock] = satpos(t, prn, nav)
% satpos: calculates satellite position and clock error from ephemeris
%
% args  : t        I   time of transmission (GPST)
%         prn      I   satellite PRN number
%         nav      I   navigation data structure (with .eph_gps and .eph_bds)
%
% return: svpos    O   satellite position in ECEF (m) [x; y; z]
%         svclock  O   satellite clock error (s)

% --- Define Physical Constants ---
MU_GPS = 3.986005E14;     % WGS-84 value of gravitational constant for GPS
OMEGA_E_GPS = 7.2921151467E-5; % WGS-84 value of earth's rotation rate
MU_BDS = 3.986004418E14;  % CGCS2000 value for BDS
OMEGA_E_BDS = 7.2921150E-5; % CGCS2000 value

svpos = [0; 0; 0];
svclock = 0;

% Determine satellite system
[sys, ~] = satsys(prn); % Using the simple version here for demonstration

if sys == 1 % SYS_GPS
    eph_list = nav.eph_gps;
    mu = MU_GPS;
    omega_e = OMEGA_E_GPS;
elseif sys == 2 % SYS_BDS
    eph_list = nav.eph_bds;
    mu = MU_BDS;
    omega_e = OMEGA_E_BDS;
else
    return; % Unsupported system
end

% Find the best ephemeris for the given time t
eph = [];
min_diff = inf;
for i = 1:length(eph_list)
    if eph_list(i).prn == prn
        diff = abs(t - eph_list(i).toe);
        if diff < min_diff
            min_diff = diff;
            eph = eph_list(i);
        end
    end
end

if isempty(eph)
    % fprintf('No ephemeris found for PRN %d at time %f\n', prn, t);
    return;
end

% --- Satellite Position and Clock Calculation (Common for GPS & BDS Keplerian part) ---
tk = t - eph.toe;
% Account for week rollover
if tk > 302400, tk = tk - 604800; end
if tk < -302400, tk = tk + 604800; end

% Mean anomaly
a = eph.sqrta^2;
n0 = sqrt(mu / a^3);
n = n0 + eph.dn;
m = eph.m0 + n * tk;

% Eccentric anomaly (iterative solution)
e = eph.e;
ek = m;
for i=1:10
    ek_old = ek;
    ek = m + e * sin(ek);
    if abs(ek - ek_old) < 1e-12, break; end
end

% True anomaly
vk = atan2(sqrt(1 - e^2) * sin(ek), cos(ek) - e);

% Argument of latitude
phik = vk + eph.omega;

% Second-order corrections
du = eph.cus * sin(2*phik) + eph.cuc * cos(2*phik);
dr = eph.crs * sin(2*phik) + eph.crc * cos(2*phik);
di = eph.cis * sin(2*phik) + eph.cic * cos(2*phik);
uk = phik + du;
rk = a * (1 - e*cos(ek)) + dr;
ik = eph.i0 + di + eph.idot * tk;

% Position in orbital plane
x_orb = rk * cos(uk);
y_orb = rk * sin(uk);

% ECEF coordinates
omega_k = eph.omega0 + (eph.omegadot - omega_e) * tk - omega_e * eph.toe;
X = x_orb * cos(omega_k) - y_orb * cos(ik) * sin(omega_k);
Y = x_orb * sin(omega_k) + y_orb * cos(ik) * cos(omega_k);
Z = y_orb * sin(ik);

svpos = [X; Y; Z];

% --- Satellite Clock Correction ---
dt = t - eph.toc;
if dt > 302400, dt = dt - 604800; end
if dt < -302400, dt = dt + 604800; end

svclock = eph.af0 + eph.af1 * dt + eph.af2 * dt^2;

% Relativistic correction
F = -2 * sqrt(mu) / (299792458^2);
svclock = svclock + F * e * eph.sqrta * sin(ek);

% --- Special BDS Corrections for GEO satellites ---
% NOTE: A complete BDS implementation would require handling GEO satellite
%       specific corrections which involve a different model.
%       This implementation uses the Keplerian model for all BDS sats,
%       which is accurate for MEO/IGSO but less so for GEO.
%       For RTK, this is often sufficient as errors are differenced out.