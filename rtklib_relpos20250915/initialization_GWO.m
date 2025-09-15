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
% �˺������ڳ�ʼ�������������ǣ��ĳ�ʼ��Ⱥ
% �������:
% SearchAgents_no: �����������ǣ�������
% dim: ������ά��
% ub: �������Ͻ磬�����ǵ�����ֵ������
% lb: �������½磬�����ǵ�����ֵ������
% �������:
% Positions: ��ʼ�������������λ�þ��󣬴�СΪ SearchAgents_no x dim
function Positions=initialization(SearchAgents_no,dim,ub,lb)

% ����߽�����������Ͻ�����������
Boundary_no= size(ub,2); 

% If the boundaries of all variables are equal and user enter a single
% number for both ub and lb
% ������б����ı߽���ͬ���û�Ϊ�Ͻ���½綼�����˵�����ֵ
if Boundary_no==1
    % ʹ�� rand ��������һ����СΪ SearchAgents_no x dim ���������
    % �������Χ�� [0, 1] ֮��
    % ͨ�� (ub - lb) ������Χ�������� lb���������ӳ�䵽 [lb, ub] ����
    Positions=rand(SearchAgents_no,dim).*(ub-lb)+lb;
end

% If each variable has a different lb and ub
% ���ÿ���������в�ͬ���Ͻ���½�
if Boundary_no>1
    % ����ÿ������ά��
    for i=1:dim
        % ��ȡ��ǰά�ȵ��Ͻ�
        ub_i=ub(i);
        % ��ȡ��ǰά�ȵ��½�
        lb_i=lb(i);
        % Ϊ��ǰά������һ����СΪ SearchAgents_no x 1 ���������
        % �������Χ�� [0, 1] ֮��
        % ͨ�� (ub_i - lb_i) ������Χ�������� lb_i���������ӳ�䵽 [lb_i, ub_i] ����
        % �����丳ֵ�� Positions ����ĵ� i ��
        Positions(:,i)=rand(SearchAgents_no,1).*(ub_i-lb_i)+lb_i;
    end
end