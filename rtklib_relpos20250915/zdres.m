function [y,e,azel]=zdres(time,obs,n,RS,DTS,rr)
% *************************************************************************
% --- VERSION 2.0: UPGRADED for new data structures and matrix dimensions ---
% *************************************************************************

global CLIGHT FREQ1 OMGE RE_WGS84; % Removed unused globals, added OMGE

% --- Initialize output arrays ---
y = zeros(2*n, 1);    % Pre-allocate for speed
e = zeros(n, 3);      % Pre-allocate for speed
azel = zeros(2*n, 1); % Pre-allocate for speed

% --- Extract observation data from the new 'obs.data' matrix ---
% The columns in obs.data are assumed to be:
% [1:PRN, 2:C1/P1, 3:L1, 4:D1, 5:S1, 6:P2, 7:L2, 8:Sys]
epoch_data = obs.data;
L1_obs = epoch_data(:, 3); % Column 3 for L1 carrier phase
C1_obs = epoch_data(:, 2); % Column 2 for C1/P1 pseudorange

rr_ = rr(1:3); % Take only the first 3 elements
rr_ = rr_(:)';   % Force it to be a 1x3 row vector
pos = ecef2pos(rr_);

for i=1:n
    % --- MODIFIED START: Correctly index the 3xn 'RS' matrix ---
    sat_pos = RS(:, i)'; % Get the i-th column and transpose to a 1x3 row vector
    
    %求星地距离和星地矢量
    e_i = sat_pos - rr_; % Vector subtraction
    r = norm(e_i);
    e(i, :) = e_i / r; % Store the unit vector
    
    %地球自转修正 (Sagnac effect)
    r = r + OMGE * (sat_pos(1)*rr_(2) - sat_pos(2)*rr_(1)) / CLIGHT;
    % --- MODIFIED END ---
    
    if(r<=0), continue; end
    
    %求仰角方位角 
    [az, el] = satelev(pos, e(i, :));
    azel(2*i-1) = az;
    azel(i*2)   = el;
    
    % Satellite clock-bias修正星地距离
    r_corrected = r - CLIGHT*DTS(i);
    
    % 对流层延迟模型修正 (Your original logic is preserved)
    dtrp = tropocorr(time, pos, el);
    r_corrected = r_corrected + dtrp;
    
    % 天线相位中心修正 (Your code sets this to zero, so we honor that)
    dant = [0,0,0];
    
    % --- MODIFIED START: Use extracted data from 'epoch_data' matrix ---
    lam1 = CLIGHT/FREQ1; % L1 wavelength
    y(2*i-1) = L1_obs(i)*lam1 - r_corrected - dant(1); % Carrier phase residual
    y(2*i)   = C1_obs(i) - r_corrected - dant(1); % Pseudorange residual
    % --- MODIFIED END ---
end

end % End of main function zdres

% --- HELPER FUNCTIONS (Refactored from your original code for clarity) ---

function [az, el] = satelev(pos, e_vec)
    global RE_WGS84;
    az = 0; el = pi/2;
    if pos(3) > -RE_WGS84
        enu = ecef2enu(pos, e_vec);
        if norm(enu(1:2)) < 1E-12
            az=0;
        else
            az=atan2(enu(1),enu(2)); 
        end
        if(az<0), az=az+2*pi; end
        el=asin(enu(3));
    end
end

function dtrp = tropocorr(time, pos, el)
    % This function encapsulates your detailed troposphere model.
    global R2D;
    coef=[1.2769934E-3, 1.2683230E-3, 1.2465397E-3, 1.2196049E-3, 1.2045996E-3;
          2.9153695E-3, 2.9152299E-3, 2.9288445E-3, 2.9022565E-3, 2.9024912E-3;
          62.610505E-3, 62.837393E-3, 63.721774E-3, 63.824265E-3, 64.258455E-3;
          0.0000000E-0, 1.2709626E-5, 2.6523662E-5, 3.4000452E-5, 4.1202191E-5;
          0.0000000E-0, 2.1414979E-5, 3.0160779E-5, 7.2562722E-5, 11.723375E-5;
          0.0000000E-0, 9.0128400E-5, 4.3497037E-5, 84.795348E-5, 170.37206E-5;
          5.8021897E-4, 5.6794847E-4, 5.8118019E-4, 5.9727542E-4, 6.1641693E-4;
          1.4275268E-3, 1.5138625E-3, 1.4572752E-3, 1.5007428E-3, 1.7599082E-3;
          4.3472961E-2, 4.6729510E-2, 4.3908931E-2, 4.4626982E-2, 5.4736038E-2];
    aht=[2.53E-5, 5.49E-3, 1.14E-3];

    zhd = 0;
    if ~(pos(3)<-100 || pos(3)>1E+4 || el<=0)
        hgt = max(0, pos(3));
        pres=1013.25*(1.0-2.2557E-5*hgt)^5.2568;
        temp=15.0-6.5E-3*hgt+273.16;
        Edtrp=6.108*0*exp((17.15*temp-4684.0)/(temp-38.45)); % Note: your original code has a '0' here, maybe a placeholder?
        z=pi/2.0-el;
        trph=0.0022768*pres/(1.0-0.00266*cos(2.0*pos(1))-0.00028*hgt/1E+3);
        if cos(z) > 1e-6, trph = trph / cos(z); end % Avoid division by zero
        trpw=0.002277*(1255.0/temp+0.05)*Edtrp;
        if cos(z) > 1e-6, trpw = trpw / cos(z); end
        zhd=trph+trpw;
    end
    
    tropmapF=0;
    if ~(pos(3)<-1000.0 || pos(3)>20000.0 || el <= 0)
        lat=pos(1)*R2D;
        y_t=(time2doy(time)-28.0)/365.25 + (1-sign(lat))/2;
        cosy=cos(2*pi*y_t);
        lat=abs(lat);
        ah=zeros(3,1); aw=zeros(3,1);
        for w=1:3
            interPC1=interpc(coef(w,:),lat);
            interPC2=interpc(coef(w+3,:),lat);
            interPC3=interpc(coef(w+6,:),lat);
            ah(w)=interPC1-interPC2*cosy;
            aw(w)=interPC3;
        end
        dm=(1/sin(el)-mapf(el,aht(1),aht(2),aht(3)))*pos(3)/1000;
        tropmapF=mapf(el,ah(1),ah(2),ah(3))+dm;
    end
    dtrp = tropmapF*zhd;
end

function interPC=interpc(coef_t,lat)
    i_t=fix(lat/15);
    if i_t<1, interPC=coef_t(1);
    elseif i_t>4, interPC=coef_t(5);
    else
        interPC=coef_t(i_t)*(1-lat/15+i_t)+coef_t(i_t+1)*(lat/15.0-i_t);
    end
end

function mapF=mapf(el_t,ah,bh,ch)
    sinel=sin(el_t);
    mapF=(1+ah/(1+bh/(1+ch)))/(sinel+(ah/(sinel+bh/(sinel+ch))));
end

function doy = time2doy(time)
    % A simple helper to get day of year from gpst
    % This might need a more robust implementation depending on your time system
    doy = floor(time / 86400) + 1; % Simple conversion, might be off by one day
end