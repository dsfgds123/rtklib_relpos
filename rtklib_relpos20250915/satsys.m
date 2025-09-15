function [sys, prn_out] = satsys(prn, sys_in)
% satsys: determine satellite system from PRN number
%
% args  : prn      I   PRN number (for RINEX 2) or satellite number (for RINEX 3)
%         sys_in   I   (optional) system character from RINEX 3 (e.g., 'G', 'C')
%
% return: sys      O   satellite system (SYS_GPS, SYS_BDS, etc.)
%         prn_out  O   PRN number within its constellation (e.g. C02 -> 2)
%
% notes : This function defines the system constants used throughout the project.

% --- Define System Constants ---
SYS_NONE = 0;
SYS_GPS  = 1;
SYS_BDS  = 2;
% Future expansion:
% SYS_GAL  = 3; % Galileo
% SYS_GLO  = 4; % GLONASS

% --- Define PRN ranges ---
MIN_PRN_GPS = 1; MAX_PRN_GPS = 32;
MIN_PRN_BDS = 1; MAX_PRN_BDS = 37; % Covering BDS-2 and BDS-3

sys = SYS_NONE;

if nargin < 2 || isempty(sys_in)
    % Fallback for RINEX 2.x or when system character is not provided
    if prn >= MIN_PRN_GPS && prn <= MAX_PRN_GPS
        sys = SYS_GPS;
    end
    prn_out = prn;
else
    % Logic for RINEX 3.x with system character
    prn_out = str2double(prn);
    switch upper(sys_in)
        case 'G'
            if prn_out >= MIN_PRN_GPS && prn_out <= MAX_PRN_GPS
                sys = SYS_GPS;
            end
        case 'C'
            if prn_out >= MIN_PRN_BDS && prn_out <= MAX_PRN_BDS
                sys = SYS_BDS;
            end
        % Add cases for 'E', 'R' etc. for future expansion
    end
end