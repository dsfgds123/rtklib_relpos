function [ns,sat,ir,ib]=selsat(obsr,nr,obsu,nb,azel_b)
% *************************************************************************
% --- VERSION 2.0: UPGRADED for new observation data structures ---
% *************************************************************************
% input
% obsr :                移动站观测数据结构体 (NEW structure)
% nr :                  移动站卫星数
% obsu :                基站观测数据结构体 (NEW structure)
% nb :                  基站卫星数
% azel_b ：             基站仰角方位角 (azel_b(2*j) is elevation for j-th base sat)
% **************************
% output
% ns ：                 共视卫星数            
% sat ：                共视卫星号 (PRNs of common satellites)
% ir ：                 移动站共视卫星在其列表中的索引
% ib ：                 基站共视卫星在其列表中的索引
% *************************

% --- Initialize output arrays ---
sat = [];
ir = [];
ib = [];
ns = 0;

% --- MODIFIED START: Extract PRN lists from the new data structures ---
% In the new structure, satellite PRNs are in the first column of the .data matrix
rover_prns = obsr.data(:, 1);
base_prns  = obsu.data(:, 1);
% --- MODIFIED END ---

% --- Your original nested loop logic is preserved, but adapted ---
for i=1:nr
    for j=1:nb
        % Compare PRNs from the extracted lists
        if rover_prns(i) == base_prns(j)
            
            % Check base station satellite elevation threshold (your logic is unchanged)
            % 0.261799... rad is approximately 15 degrees.
            if azel_b(2*j) >= 0.34906585%0.261799387799149 %15°
                ns = ns + 1;
                sat(ns) = rover_prns(i); % Store the common PRN
                ir(ns) = i;              % Store the index from the rover list
                ib(ns) = j;              % Store the index from the base list
                break; % Found a match for this rover satellite, move to the next one
            end
            
        end
    end
end

% Note: Your original code had an issue where 'ns' could be 1 even if no
% satellites were found because it was initialized to 1. The new logic
% initializes ns=0 and only increments it when a valid common satellite is found,
% which is more robust. If your old code `ns=ns-1` was intended to fix this,
% the new logic solves the root cause.