function [xp_,Pp_]=EKF_filter(x,P,H,v,R,n,m,HOLD)
% *************************************************************************
% --- VERSION 2.0: Hardened with robustness checks and corrected vector orientation ---
% *************************************************************************
% input
% x :   State vector (nx1 column vector)
% P :   State covariance matrix (nxn)
% H :   Design matrix (mxn)
% v :   Innovation/residual vector (mx1 column vector)
% R :   Measurement covariance matrix (mxm)
% n :   Number of states (rtk.nx)
% m :   Number of measurements (nv)
% HOLD: (not used in this logic, but kept for compatibility)
% ***************************
% output
% xp_ : Updated state vector (nx1 column vector)
% Pp_ : Updated state covariance matrix (nxn)
% ***************************

% Initialize output with input values. If update fails, we return the original state.
xp_=x;
Pp_=P;

% --- Find indices of active states ---
ix = zeros(n, 1); % Pre-allocate
k=0;
for i=1:n
    if x(i)~=0 && P(i,i)>0 % A more standard check for active states
        k=k+1;
        ix(k)=i;
    end
end
if k > 0
    ix = ix(1:k); % Trim to actual size
else
    % --- MODIFIED START: Robustness check for zero active states ---
    % If k=0, it means no states are active for update. Do nothing and return.
    return;
    % --- MODIFIED END ---
end

% --- Create subset matrices for active states ---
x_subset = x(ix); % Creates a kx1 column vector
P_subset = P(ix, ix); % Creates a kxk covariance matrix
H_subset = H(:, ix); % Creates an mxk design matrix

% --- Kalman Filter Update Equations ---
% Note: The state 'x_subset' is already a column vector.
% The innovation 'v' should also be a column vector.
% The original 'v=v'' and 'x_=x_'' are removed as they caused dimension issues.

% Kalman Gain K = P*H' * inv(H*P*H' + R)
K = P_subset * H_subset' / (H_subset * P_subset * H_subset' + R);

% State Update: x_new = x_old + K * innovation
x_updated = x_subset + K * v;

% Covariance Update: P_new = (I - K*H) * P_old
I = eye(k);
P_updated = (I - K * H_subset) * P_subset;

% --- Re-embed the updated subset back into the full state vector and covariance matrix ---
xp_(ix) = x_updated;
Pp_(ix, ix) = P_updated;

end