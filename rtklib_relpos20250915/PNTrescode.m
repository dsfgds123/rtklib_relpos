function [v,var,nv,H]=PNTrescode(X,t,rs,dts,vare,Snum,R,Code,nav,Sys)
% *************************************************************************
% --- CODE UPGRADED to handle 3xSnum 'rs' matrix and new 'nav' structure ---
% *************************************************************************

% --- MODIFIED START: Added OMGE to the global declaration ---
global CLIGHT FREQ1 OMGE;
% --- MODIFIED END ---

H=[]; 
v = [];
var = [];
rr = X(1:3);
dtr = X(4);
nv=0;

pos=ecef2pos(rr);

for j=1:Snum
    sat_pos = rs(:, j);
    
    e = sat_pos - rr;
    r = norm(e);
    
    % Sagnac effect correction (This line will now work correctly)
    r = r + OMGE * (sat_pos(1)*rr(2) - sat_pos(2)*rr(1)) / CLIGHT;
    
    if(r<=0), continue; end
    
    e = e / r;
    
    [az, el] = satelev(pos, e);
    
    if(el < (15*pi/180.0)), continue; end
    
    prn_j = Code(j);
    sys_j = Sys(j);
    tgd = 0;
    
    if sys_j == 1 % GPS
        eph_list = nav.eph_gps;
    elseif sys_j == 2 % BDS
        eph_list = nav.eph_bds;
    else
        eph_list = [];
    end

    if ~isempty(eph_list)
        for eph_idx = 1:length(eph_list)
            if eph_list(eph_idx).prn == prn_j
                if isfield(eph_list(eph_idx), 'tgd') && ~isempty(eph_list(eph_idx).tgd)
                    tgd = CLIGHT * eph_list(eph_idx).tgd(1);
                else
                    tgd = 0; % Ensure tgd is defined if the field doesn't exist
                end
                break;
            end
        end
    end
    
    % In RINEX 3, TGD is usually for GPS L1/L2. B1/B2 TGD is often TGD1.
    % For GPS+BDS SPP, a common simplification is to only apply GPS TGD and ignore BDS TGD (BDS D1/D2)
    % as its effect is small and often absorbed by the receiver clock estimate.
    if sys_j == 1
        P1_P2 = tgd; 
    else
        P1_P2 = 0;
    end
    PC = R(j) - P1_P2;
    vmeas = (0.3)^2;
    
    [dion, vion] = ionocorr(t, pos, [az, el]);
    [dtrp, vtrp] = tropocorr(pos, el);
    
    current_v = PC - (r + dtr - CLIGHT*dts(j) + dion + dtrp);
    current_H = [-e(1), -e(2), -e(3), 1];
    current_var = (100^2)*((0.003^2)+(0.003^2)/sin(el)) + vare(j) + vmeas + vion + vtrp;
    
    v = [v; current_v];
    var = [var; current_var];
    H = [H; current_H];
    nv = nv + 1;
end

end

% --- Helper functions (Please ensure these are populated with your original code) ---
function [az, el] = satelev(pos, e)
    RE_WGS84 = 6378137.0;
    if pos(3) > -RE_WGS84
        enu = ecef2enu(pos, e);
        if norm(enu(1:2)) < 1e-12, az = 0; else, az = atan2(enu(1), enu(2)); end
        if az < 0, az = az + 2*pi; end
        el = asin(enu(3));
    else
        az = 0; el = pi/2;
    end
end

function [dion, vion] = ionocorr(t, pos, azel)
    % (Paste your original ionosphere correction code here)
    dion = 0; vion = (5.0)^2;
end

function [dtrp, vtrp] = tropocorr(pos, el)
    % (Paste your original troposphere correction code here)
    dtrp = 0; vtrp = (0.3)^2;
end