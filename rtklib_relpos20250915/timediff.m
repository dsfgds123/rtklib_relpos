function t = timediff(t1, t2)
% timediff: computes the difference between two time epochs
% --- VERSION 2.0: UPGRADED to handle scalar GPST values ---
%
% args  : t1, t2   I   Time epochs as scalar GPST (seconds of week)
% return: t        O   Time difference in seconds (t1 - t2)

% The new time format is a simple scalar double, so the difference
% is a simple subtraction. This is more robust and standard in MATLAB.
t = t1 - t2;

% The original code was: t=(t1(1)-t2(1))+t1(2)-t2(2);
% This was designed for an obsolete 2-element time format (e.g., [integer_sec, frac_sec])
% and is no longer needed.

end