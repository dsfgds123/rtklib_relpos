function sol=pointpos(sol,obs,nav)
% *************************************************************************
% --- FINAL VERSION: Updates the input 'sol' struct, preserving its fields ---
% *************************************************************************

% Extract necessary info from input observation
tr = obs.time_vec;
epoch_data = obs.data;
Snum = size(epoch_data, 1);

if Snum < 4 
    sol.stat = 0; % SOLQ_NONE
    return;
end

Code = epoch_data(:, 1);
R = epoch_data(:, 2);
Sys = epoch_data(:, 8);

tJD=TimetoJD(tr(1),tr(2),tr(3),tr(4),tr(5),tr(6));

[rs,dts,vare] = satposs(obs.time, tJD, Code, Snum, nav, R, Sys);

% --- MODIFIED: Call the new estpos and update the sol struct ---
% estpos now returns raw values, not a struct
[rr_new, dtr_new, stat_new] = estpos(obs.time, sol.rr, rs, dts, Snum, R, Code, vare, nav, Sys);

% Update the fields of the existing sol struct
sol.time = obs.time;
sol.rr = rr_new;
sol.dtr = dtr_new;
sol.stat = stat_new;
% Note: other fields like .qr, .ratio, .ns are untouched, preserving the structure
% --- MODIFIED END ---

end