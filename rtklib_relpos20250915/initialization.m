%___________________________________________________________________%
%  Grey Wolf Optimizer (GWO) source codes version 1.0               %
%                                                                   %
%  Developed in MATLAB R2011b(7.13)                                 %
%                                                                   %
%  Author and programmer: Seyedali Mirjalili                        %
%                                                                   %
%         e-Mail: ali.mirjalili@gmail.com                           %
%                 seyedali.mirjalili@griffithuni.edu.au             %
%                                                                   %
%       Homepage: http://www.alimirjalili.com                       %
%                                                                   %
%   Main paper: S. Mirjalili, S. M. Mirjalili, A. Lewis             %
%               Grey Wolf Optimizer, Advances in Engineering        %
%               Software , in press,                                %
%               DOI: 10.1016/j.advengsoft.2013.12.007               %
%                                                                   %
%___________________________________________________________________%

% This function initialize the first population of search agents
% 此函数用于初始化搜索代理（灰狼）的初始种群
% 输入参数:
% SearchAgents_no: 搜索代理（灰狼）的数量
% dim: 变量的维度
% ub: 变量的上界，可以是单个数值或向量
% lb: 变量的下界，可以是单个数值或向量
% 输出参数:
% Positions: 初始化后的搜索代理位置矩阵，大小为 SearchAgents_no x dim
function Positions=initialization(SearchAgents_no,dim,ub,lb)

% 计算边界的数量，即上界向量的列数
Boundary_no= size(ub,2); 

% If the boundaries of all variables are equal and user enter a single
% number for both ub and lb
% 如果所有变量的边界相同，用户为上界和下界都输入了单个数值
if Boundary_no==1
    % 使用 rand 函数生成一个大小为 SearchAgents_no x dim 的随机矩阵
    % 随机数范围在 [0, 1] 之间
    % 通过 (ub - lb) 调整范围，并加上 lb，将随机数映射到 [lb, ub] 区间
    Positions=rand(SearchAgents_no,dim).*(ub-lb)+lb;
end

% If each variable has a different lb and ub
% 如果每个变量都有不同的上界和下界
if Boundary_no>1
    % 遍历每个变量维度
    for i=1:dim
        % 获取当前维度的上界
        ub_i=ub(i);
        % 获取当前维度的下界
        lb_i=lb(i);
        % 为当前维度生成一个大小为 SearchAgents_no x 1 的随机向量
        % 随机数范围在 [0, 1] 之间
        % 通过 (ub_i - lb_i) 调整范围，并加上 lb_i，将随机数映射到 [lb_i, ub_i] 区间
        % 并将其赋值给 Positions 矩阵的第 i 列
        Positions(:,i)=rand(SearchAgents_no,1).*(ub_i-lb_i)+lb_i;
    end
end