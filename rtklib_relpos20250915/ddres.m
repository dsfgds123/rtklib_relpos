function [nv,rtk,v,H,R,vflg]=ddres(rtk,dt,x,sat,y_r,e_r,azel_r,y_b,ir,ib,ns,Hmode)
% *************************************************************************
% --- VERSION 2.0: Hardened with robustness checks ---
% *************************************************************************

global CLIGHT FREQ1;

% --- Initialize outputs to safe default values ---
nv = 0;
v = [];
H = [];
R = [];
vflg = [];

% --- Robustness Check: Ensure there are enough satellites to form a double difference ---
if ns < 2
    return; % Not enough common satellites to proceed
end

% Pre-allocate matrices for speed. A double difference pair is formed for each non-ref sat.
% This is done for phase (f=1) and code (f=2). So max nv is 2*(ns-1).
max_nv = 2 * (ns - 1);
v = zeros(max_nv, 1);
H = zeros(max_nv, rtk.nx);
vflg = zeros(max_nv, 1);
Ri = zeros(max_nv, 1);
Rj = zeros(max_nv, 1);

% Your original logic for baseline calculation
[bl,dr]=baseline(x,rtk.rb);

% --- Reset residuals in satellite status ---
for i=1:length(sat)
    rtk.ssat(sat(i)).resp = 0;
    rtk.ssat(sat(i)).resc = 0;
end

% Loop for phase (f=1) and code (f=2) observations
nb = [0, 0];
for f=1:2
    % --- Search for reference satellite with the highest elevation ---
    ref_idx = -1; % Index within the 'ns' list (1 to ns)
    max_el = 0;
    for j=1:ns
        el = azel_r(2*ir(j)); % Elevation of j-th common satellite
        if el > max_el
            max_el = el;
            ref_idx = j;
        end
    end
    
    % --- Robustness Check: Ensure a valid reference satellite was found ---
    if ref_idx < 1
        continue; % Skip this frequency if no valid reference satellite
    end
    
    % The index 'i' now represents the reference satellite
    i = ref_idx;
    
    % --- Make double differences against the reference satellite ---
    for j=1:ns
        if i == j, continue; end % Skip differencing with itself

        lami = CLIGHT/FREQ1;
        lamj = CLIGHT/FREQ1;
        
        % Double-differenced residual (your formula is correct and preserved)
        current_v = (y_r(f + (ir(i)-1)*2) - y_b(f + (ib(i)-1)*2)) ...
                  - (y_r(f + (ir(j)-1)*2) - y_b(f + (ib(j)-1)*2));
        
        % Partial derivatives by rover position
        if Hmode==1
            current_H_row = zeros(1, rtk.nx);
            current_H_row(1:3) = -e_r(ir(i),:) + e_r(ir(j),:);
        end
        
        % Double-differenced phase-bias term
        if f==1
            % Correct the residual for ambiguity terms
            current_v = current_v - lami*rtk.x(3+sat(i)) + lamj*rtk.x(3+sat(j));
            if Hmode==1
                % Update the H matrix for ambiguity terms
                current_H_row(3+sat(i)) = lami;
                current_H_row(3+sat(j)) = -lamj;
            end
        end
        
        % Store residuals in the rtk structure
        if f==1
            rtk.ssat(sat(j)).resc = current_v;
        else
            rtk.ssat(sat(j)).resp = current_v;
        end
        
        % Check innovation threshold (your logic preserved)
        if abs(current_v) > 30
            if f==1
                rtk.ssat(sat(i)).rejc = rtk.ssat(sat(i)).rejc + 1;
                rtk.ssat(sat(j)).rejc = rtk.ssat(sat(j)).rejc + 1;
            end
             % This was commented out in your code, so I keep it that way.
             % If uncommented, it might cause the 'Ri not found' error.
             % continue; 
        end
        
        % If we've reached here, the data is valid for this pair.
        nv = nv + 1;

        % Store the calculated values into pre-allocated arrays
        v(nv) = current_v;
        if Hmode==1
            H(nv,:) = current_H_row;
        end
        
        % Single-differenced measurement error variances (your vareer function is used)
        Ri(nv) = vareer(azel_r(ir(i)*2), bl, dt, f);
        Rj(nv) = vareer(azel_r(ir(j)*2), bl, dt, f);
        
        % Set valid data flags
        if f==1
            rtk.ssat(sat(i)).vsat = 1;
            rtk.ssat(sat(j)).vsat = 1;
        end
        
        % vflg packing (your logic preserved)
        obs_type_flag = (f==2); % 0 for phase, 1 for code
        vflg(nv) = bitor(bitor(bitor(bitshift(sat(i),16), bitshift(sat(j),8)), bitshift(obs_type_flag,4)), 0);
        
        nb(f) = nb(f) + 1;
    end
end

% Trim unused parts of the pre-allocated arrays
if nv > 0
    v = v(1:nv);
    H = H(1:nv, :);
    vflg = vflg(1:nv);
    Ri = Ri(1:nv);
    Rj = Rj(1:nv);
    
    % Double-differenced measurement error covariance (your ddcov is used)
    R = ddcov(nb, f, Ri, Rj, nv);
else
    % If no valid double differences were formed, return empty matrices
    v = []; H = []; R = []; vflg = [];
end

end % End of main function ddres

% --- Your original helper functions, unchanged ---
function [bl_,dr_]=baseline(rr,rb)
    dr_ = rr - rb;
    bl_ = norm(dr_);
end

function Ri_nv=vareer(el,bl,dt,f)
    CLIGHT=299792458.0;
    c=0*bl/1e4;
    d=CLIGHT*5e-12*dt;
    fact = 1;
    if f==2, fact=100; end
    a=fact*0.003;
    b=fact*0.003;
    Ri_nv=2*(a*a+b*b/sin(el)/sin(el)+c*c)+d*d;
end

function R_=ddcov(nb_,b_,Ri_,Rj_,nv_)
    R_=zeros(nv_);
    k=0;
    for b=1:b_ % b_ is the number of frequencies (1 or 2)
        n_freq = nb_(b); % Number of DD pairs for this frequency
        for i=1:n_freq
            for j=1:n_freq
                if i==j
                    R_(k+i,k+j) = Ri_(k+i) + Rj_(k+i);
                else
                    % Covariance term is the variance of the common reference satellite
                    R_(k+i,k+j) = Ri_(k+i); % Or Ri_(k+j), they should be the same
                end
            end
        end
        k = k + n_freq;
    end
end