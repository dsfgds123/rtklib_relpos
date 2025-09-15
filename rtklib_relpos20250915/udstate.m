function rtk=udstate(mode,ARmode,rtk,obsr,obsu,sat,ir,ib,ns,nav)
% *************************************************************************
% --- VERSION 2.0: Adapted for new data structures and multi-constellation ---
% *************************************************************************

global MAXSAT;

VAR_POS = 900.0;
tt = abs(rtk.tt);

rtkX = rtk.x;
rtkP = rtk.P;

% --- (1) Update position state ---
% Initialize position for the first epoch
if norm(rtkX(1:3)) <= 0
    for i=1:3
        [rtkX,rtkP] = initx(rtkX, rtkP, rtk.sol.rr(i), VAR_POS, i);
    end
else
    % For kinematic mode, re-initialize position at every epoch
    % For static mode, position propagates from the previous epoch
    if mode == 0 % Kinematic
        for i=1:3
            [rtkX,rtkP] = initx(rtkX, rtkP, rtk.sol.rr(i), VAR_POS, i);
        end
    end
end

% --- (2) Update phase-bias state ---

% Detect cycle slips for currently visible satellites
for i=1:ns
    sat_prn = sat(i);
    % Your original logic for cycle slip detection seems to rely on LLI flags
    % which we haven't implemented in readrinexobs3. For now, this will be inactive.
    % To make this work, 'detslp_ll' would need to be updated.
    % rtk.ssat(sat_prn).slip = detslp_ll(rtk, obsr, ir(i), 1);
    % rtk.ssat(sat_prn).slip = detslp_ll(rtk, obsu, ib(i), 2);
end

% Reset phase-bias for all possible satellites if conditions are met
for i=1:MAXSAT
    rtk.ssat(i).outc = rtk.ssat(i).outc + 1;
    reset = (rtk.ssat(i).outc > 5); % 5 is outage threshold
    
    amb_idx = 3 + i;
    if ARmode == 0 && rtkX(amb_idx) ~= 0
        [rtkX,rtkP] = initx(rtkX, rtkP, 0, 0, amb_idx);
    elseif reset && rtkX(amb_idx) ~= 0
        [rtkX,rtkP] = initx(rtkX, rtkP, 0, 0, amb_idx);
    end
    
    if ARmode ~= 0 && reset
        rtk.ssat(i).lock = 0;
    end
end

% Correct phase-bias offset
bias = zeros(ns,1);
offset = 0;
j_offset = 0;

for i=1:ns
    sat_prn = sat(i);
    [cp, pr] = sgobs(obsr, obsu, ir(i), ib(i)); % Get phase and code diff
    
    if cp==0 || pr==0, continue; end
    
    amb_idx = 3 + sat_prn;
    bias(i) = cp - pr / (CLIGHT/FREQ1);
    
    if rtkX(amb_idx) ~= 0
        offset = offset + bias(i) - rtkX(amb_idx);
        j_offset = j_offset + 1;
    end
end

% Apply offset correction if applicable
if j_offset > 0
    offset_mean = offset / j_offset;
    for i=1:MAXSAT
        amb_idx = 3 + i;
        if rtkX(amb_idx) ~= 0
            rtkX(amb_idx) = rtkX(amb_idx) + offset_mean;
        end
    end
end

% Set initial states of phase-bias for new satellites
for i=1:ns
    sat_prn = sat(i);
    amb_idx = 3 + sat_prn;
    if bias(i) ~= 0 && rtkX(amb_idx) == 0
        [rtkX,rtkP] = initx(rtkX, rtkP, bias(i), 900, amb_idx);
    end
end

rtk.x = rtkX;
rtk.P = rtkP;

end % End of main function udstate

% --- Helper function initx (assuming it takes these inputs) ---
function [x,P] = initx(x, P, val, var, index)
    x(index) = val;
    for i=1:length(P)
        P(i, index) = 0;
        P(index, i) = 0;
    end
    P(index, index) = var;
end

% --- Helper function sgobs, updated for the new data structure ---
function [cp,pr] = sgobs(obsr, obsu, ir, ib)
    global CLIGHT FREQ1;

    % Columns: [1:PRN, 2:C1/P1, 3:L1, ...]
    L1_r = obsr.data(ir, 3);
    L1_b = obsu.data(ib, 3);
    C1_r = obsr.data(ir, 2);
    C1_b = obsu.data(ib, 2);

    if L1_r==0 || L1_b==0, cp=0; else, cp = L1_r - L1_b; end
    if C1_r==0 || C1_b==0, pr=0; else, pr = C1_r - C1_b; end
end