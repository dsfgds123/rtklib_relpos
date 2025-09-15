function [rs, dts, vare] = satposs(time, tJD, Code, Snum, nav, R, Sys)
% satposs: calculates satellite positions, clock errors, and variances from ephemeris
% --- VERSION 2.0: UPGRADED FOR GPS+BDS and new input arguments ---
%
% args  : time     I   time of transmission (GPST seconds of week)
%         tJD      I   Julian Day for time of transmission
%         Code     I   Vector of satellite PRNs for this epoch
%         Snum     I   Number of satellites in this epoch
%         nav      I   Navigation data structure (with .eph_gps and .eph_bds)
%         R        I   Vector of pseudoranges (not directly used here, but kept for compatibility)
%         Sys      I   Vector of satellite system IDs (1=GPS, 2=BDS)
%
% return: rs       O   Satellite positions in ECEF (m) [3 x Snum]
%         dts      O   Satellite clock errors (s) [1 x Snum]
%         vare     O   Satellite variances (m^2) [1 x Snum]

% --- Define Physical Constants ---
MU_GPS = 3.986005E14;     % WGS-84 value of gravitational constant for GPS
OMEGA_E_GPS = 7.2921151467E-5; % WGS-84 value of earth's rotation rate
MU_BDS = 3.986004418E14;  % CGCS2000 value for BDS
OMEGA_E_BDS = 7.2921150E-5; % CGCS2000 value
C = 299792458; % Speed of light

% --- Initialize Output Arrays ---
rs = zeros(3, Snum);
dts = zeros(1, Snum);
vare = zeros(1, Snum);

% --- Loop through each satellite observed in this epoch ---
for i = 1:Snum
    prn = Code(i);
    sys = Sys(i);
    
    eph_list = [];
    mu = 0;
    omega_e = 0;
    
    % Select the correct ephemeris list and constants based on system
    if sys == 1 % SYS_GPS
        if isfield(nav, 'eph_gps') && ~isempty(nav.eph_gps)
            eph_list = nav.eph_gps;
            mu = MU_GPS;
            omega_e = OMEGA_E_GPS;
        end
    elseif sys == 2 % SYS_BDS
        if isfield(nav, 'eph_bds') && ~isempty(nav.eph_bds)
            eph_list = nav.eph_bds;
            mu = MU_BDS;
            omega_e = OMEGA_E_BDS;
        end
    end
    
    if isempty(eph_list)
        continue; % Skip if no ephemeris for this system
    end

    % Find the best ephemeris for the given time
    eph = [];
    min_diff = inf;
    for j = 1:length(eph_list)
        if eph_list(j).prn == prn
            diff = abs(time - eph_list(j).toe);
            if diff < min_diff
                min_diff = diff;
                eph = eph_list(j);
            end
        end
    end

    if isempty(eph)
        continue; % Skip if no ephemeris found for this satellite
    end
    
    % --- Satellite Position and Clock Calculation (Common for GPS & BDS) ---
    tk = time - eph.toe;
    % Account for week rollover
    if tk > 302400, tk = tk - 604800; end
    if tk < -302400, tk = tk + 604800; end

    a = eph.sqrta^2;
    n0 = sqrt(mu / a^3);
    n = n0 + eph.dn;
    m = eph.m0 + n * tk;

    e_ecc = eph.e;
    ek = m;
    for k_iter=1:10
        ek_old = ek;
        ek = m + e_ecc * sin(ek);
        if abs(ek - ek_old) < 1e-12, break; end
    end

    vk = atan2(sqrt(1 - e_ecc^2) * sin(ek), cos(ek) - e_ecc);
    phik = vk + eph.omega;

    du = eph.cus * sin(2*phik) + eph.cuc * cos(2*phik);
    dr = eph.crs * sin(2*phik) + eph.crc * cos(2*phik);
    di = eph.cis * sin(2*phik) + eph.cic * cos(2*phik);
    uk = phik + du;
    rk = a * (1 - e_ecc*cos(ek)) + dr;
    ik = eph.i0 + di + eph.idot * tk;

    x_orb = rk * cos(uk);
    y_orb = rk * sin(uk);

    omega_k = eph.omega0 + (eph.omegadot - omega_e) * tk - omega_e * eph.toe;
    X = x_orb * cos(omega_k) - y_orb * cos(ik) * sin(omega_k);
    Y = x_orb * sin(omega_k) + y_orb * cos(ik) * cos(omega_k);
    Z = y_orb * sin(ik);

    rs(:, i) = [X; Y; Z];

    % --- Satellite Clock Correction ---
    dt = time - eph.toc;
    if dt > 302400, dt = dt - 604800; end
    if dt < -302400, dt = dt + 604800; end

    svclock = eph.af0 + eph.af1 * dt + eph.af2 * dt^2;
    % Relativistic correction
    F = -2 * sqrt(mu) / (C^2);
    svclock = svclock + F * e_ecc * eph.sqrta * sin(ek);
    
    dts(i) = svclock;
    
    % --- Variance Calculation (simplified) ---
    % A simple model for user range error (URE) variance
    vare(i) = 2.0^2; % Assuming a constant URE of 2 meters
end

end % End of function satposs