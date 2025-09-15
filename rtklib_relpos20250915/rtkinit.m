function rtk=rtkinit
% *************************************************************************
% --- VERSION 2.0: Dynamically sized for multi-constellation (GPS+BDS) ---
% *************************************************************************

global MAXSAT;

% --- Robustness check: Ensure MAXSAT is defined and large enough ---
if isempty(MAXSAT) || MAXSAT < 32
    fprintf('Warning: Global MAXSAT not set or too small. Defaulting to 60.\n');
    MAXSAT = 60;
end

rtk=struct;
rtk.nfix=0;

% --- Solution structure ---
rtk.sol.rr    = [0,0,0];
rtk.sol.qr    = [0,0,0];
rtk.sol.dtr   = 0;
rtk.sol.stat  = 0;
rtk.sol.ns    = 0;
rtk.sol.age   = 0;
rtk.sol.ratio = 0;
rtk.sol.time  = [0,0];

% --- Base station and time info ---
rtk.rb=[0,0,0];
rtk.tt=0;

% --- MODIFIED START: Dynamic state vector sizing ---
% State vector is [x, y, z, Amb_1, Amb_2, ..., Amb_MAXSAT]'
rtk.nx = 3 + MAXSAT; % Total number of states
rtk.na = 3; % Number of position states
% --- MODIFIED END ---

% --- Initialize state vector and covariance matrices with correct size ---
rtk.x  = zeros(rtk.nx, 1);
rtk.P  = zeros(rtk.nx, rtk.nx);
rtk.xa = zeros(rtk.na, 1);
rtk.Pa = zeros(rtk.na, rtk.na);

% --- Initialize satellite status structure for all possible satellites ---
for i=1:MAXSAT
    rtk.ssat(i).resp = 0; % Pseudorange DD residual
    rtk.ssat(i).resc = 0; % Carrier phase DD residual
    rtk.ssat(i).vsat = 0; % Valid satellite flag
    rtk.ssat(i).fix  = 0; % Fix status (0:float, 1:fix, 2:hold)
    rtk.ssat(i).slip = 0; % Cycle slip flag
    rtk.ssat(i).lock = 0; % Lock counter
    rtk.ssat(i).outc = 0; % Outage counter
    rtk.ssat(i).rejc = 0; % Rejection counter
    rtk.ssat(i).azel = [0,0];
end

end