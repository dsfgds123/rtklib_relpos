% function relpos
% ... (version history) ...
% --- CODE UPGRADED on [Current Date] to support RINEX 3.02 (GPS+BDS) ---
% --- FINAL VERSION with robust epoch synchronization and correct struct handling ---

close all;clear ;clc;
tic

% --- Global Constants ---
global CLIGHT SOLQ_NONE SOLQ_SINGLE SOLQ_FLOAT SOLQ_FIX MINFIX MAXSAT;
global FREQ1 FREQ2 D2R R2D RE_WGS84 FE_WGS84 OMGE;

CLIGHT=299792458;
SOLQ_NONE=0; SOLQ_SINGLE=1; SOLQ_FLOAT=2; SOLQ_FIX=3;
MINFIX=10;
MAXSAT=60;
FREQ1=1.57542E+9; FREQ2=1.22760E+9;
D2R=pi/180.0; R2D=180.0/pi;
RE_WGS84=6378137.0; FE_WGS84=(1.0/298.257223563);
OMGE=7.2921151467E-5;

% --- Data Loading ---
fprintf('Reading RINEX 3.02 Navigation file...\n');
nav = readrinexnav3('C:\Users\tjc\Desktop\读取gps加北斗改进经历\rtklib_relpos20250915\1h数据\base.nav'); % Use your actual RINEX 3 nav file

fprintf('Reading RINEX 3.02 Rover OBS file...\n');
ObsODat_rover = readrinexobs3('C:\Users\tjc\Desktop\读取gps加北斗改进经历\rtklib_relpos20250915\1h数据\rover.obs'); % Use your actual RINEX 3 rover obs file

fprintf('Reading RINEX 3.02 Base OBS file...\n');
ObsODat_base = readrinexobs3('C:\Users\tjc\Desktop\读取gps加北斗改进经历\rtklib_relpos20250915\1h数据\base.obs'); % Use your actual RINEX 3 base obs file

% --- Approximate Coordinates ---
Xk0=-1178093.1570+2.1130;
Yk0=5838582.6346+0.0233;
Zk0=2275048.8224+0.0242;

% --- RTK Initialization ---
rtk=rtkinit;
ARmode=2;   % 2: fix and hold
Mode=0;     % 0: kinematic
rtk.rb=[-1178093.1570, 5838582.6346, 2275048.8224]; % 基站坐标

% --- Robust Main Processing Loop with Epoch Synchronization ---
epoch_num_r = length(ObsODat_rover);
epoch_num_b = length(ObsODat_base);

num_epochs_to_process = min(epoch_num_r, epoch_num_b);
sol = repmat(rtk.sol, num_epochs_to_process, 1);
xchange = zeros(num_epochs_to_process, 1);
ychange = zeros(num_epochs_to_process, 1);
zchange = zeros(num_epochs_to_process, 1);
ratio   = zeros(num_epochs_to_process, 1);

kr = 1; kb = 1; k_out = 0;

fprintf('Starting epoch-by-epoch processing...\n');
while kr <= epoch_num_r && kb <= epoch_num_b

    time_r = ObsODat_rover(kr).time;
    time_b = ObsODat_base(kb).time;

    if abs(time_r - time_b) < 0.01
        k_out = k_out + 1;

        if k_out == 1
            pretime=[0,0];
        else  
            pretime=sol(k_out-1).time;
        end

        % --- MODIFIED: Pass the whole rtk.sol struct to pointpos for updating ---
        rtk.sol = pointpos(rtk.sol, ObsODat_rover(kr), nav);

        if pretime(1) ~= 0
            rtk.tt = timediff(rtk.sol.time, pretime);
        end

        rtk = relativepos(rtk, ObsODat_rover(kr), ObsODat_base(kb), nav, Mode, ARmode);

        sol(k_out) = rtk.sol; % This assignment will now work perfectly
        xchange(k_out) = Xk0-sol(k_out).rr(1);
        ychange(k_out) = Yk0-sol(k_out).rr(2);
        zchange(k_out) = Zk0-sol(k_out).rr(3);
        ratio(k_out) = rtk.sol.ratio;

        kr = kr + 1; kb = kb + 1;

    elseif time_r < time_b
        kr = kr + 1;
    else
        kb = kb + 1;
    end
end
fprintf('Processing finished. Total matched epochs: %d\n', k_out);

% --- Trim result arrays ---
if k_out > 0
    xchange = xchange(1:k_out);
    ychange = ychange(1:k_out);
    zchange = zchange(1:k_out);
    ratio = ratio(1:k_out);
end

% --- Plotting section ---
figure(1);
plot(xchange,'b.-'); xlabel('历元'); ylabel('change(m)'); hold on
plot(ychange,'r.-');
plot(zchange,'k.-');
legend('xchange','ychange','zchange');
title('Position Error (E-N-U components)');
grid on;

figure(2);
plot(ratio, 'g.-'); title('Ratio'); xlabel('历元'); ylabel('Ratio Value');
grid on;

toc