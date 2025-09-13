%function relpos
% version history ***********************************************************
% 2017/05/12
% Ϊ�򻯳���ʹ������ƶ�վ�ͻ�վ���ݿ�ʼ�����ʱ����ͬ����Ԫ����ͬ��ֱ��˫��
% 2017/05/24
% ������contineousģ���Ƚ���ģʽ��������static��̬����ģʽ
% 2017/05/28
% ������fix and hold ģ���Ƚ���ģʽ
% 2017/06/05
% �޸��˿������˲��е�bug��ѭ�������еĿ������˲���holdamb�����еĿ������˲��Ծ���H���󷨲�ͬ
% **************************************************************************

close all;clear ;clc;
tic

 % if 0
global HeadODat_rover;%�ƶ�վ��ͷ�ļ�
global ObsODat_rover; %�ƶ�վ�Ĺ۲�����
global HeadODat_base; %��վ��ͷ�ļ�
global ObsODat_base;  %��վ�Ĺ۲�����
global EphDat;  %��������

global CLIGHT;%����
global SOLQ_NONE;%�޽� 
global SOLQ_SINGLE;%�����
global SOLQ_FLOAT;%�����
global SOLQ_FIX;%�̶���
global MINFIX;%��С�̶���
global MAXSAT;%���ɽ��յ������� 
global FREQ1;%L1Ƶ��
global FREQ2;%L2Ƶ��
global D2R;
global R2D;
global RE_WGS84;%���򳤰뾶
global FE_WGS84;%����
global OMGE;%������ת���ٶȴ�С

CLIGHT=299792458;%����
SOLQ_NONE=0;
SOLQ_SINGLE=1;
SOLQ_FLOAT=2;
SOLQ_FIX=3;
MINFIX=10;
MAXSAT=32;
FREQ1=1.57542E+9;    %L1Ƶ��
FREQ2=1.22760E+9;    %L2Ƶ��
D2R=pi/180.0;%1������Ӧ�Ļ���
R2D=180.0/pi;%1��������Ӧ�Ķ�
RE_WGS84=6378137.0;%���뾶
FE_WGS84=(1.0/298.257223563);%����
OMGE=7.2921151467E-5;%������ת���ٶ�

EphDat=ReadGpsData;                                       %��N�ļ�
[HeadODat_rover,ObsODat_rover]=ReadObsData(0);            %���ƶ�վO�ļ�
[HeadODat_base,ObsODat_base]=ReadObsData(1);              %����վO�ļ�

 % end
%  load ('F:\��-����\matlab����\rtklib_relpos\data\text.mat');
% load ('F:\��-����\matlab����\rtklib_relpos\data\text2.mat');
N = size(EphDat,2);                                       %�õ�EphDat(i)��i�ĸ�����N�ļ��������ܸ�������o�ļ��У�����Ԫ��ͬ���������ܲ�ͬ��
epochNUM=size(ObsODat_rover,2);                           %�õ���Ԫ�������۲����ݵĸ�����ÿ��Ԫ�õ�һ�����ݣ�

% Xk0=-2005024.6830;                                        %��վ��Ľ�������
% Yk0=5411166.0497;
% Zk0=2707856.7714;

% Xk0=-2005009.6714;
% Yk0=5411172.1014;
% Zk0=2707856.3983; %data����λ��

Xk0=-2005033.7825+2.1130;%enu
Yk0=5411163.2342+0.0233;
Zk0=2707856.6419+0.0242;
%��������Ϊdata2�ļ����е�rinex2.10�ļ���
% Xk0=-2005035.76699731;
% Yk0=5411162.50934514;
% Zk0=2707856.67218072;
%��������Ϊdata2�ļ����е�a.nav��base��rover
% Xk0=-1991315.213867; 
% Yk0=5418624.208629; 
% Zk0=2702815.522514;



rtk=rtkinit;%��ʼ��
ARmode=2;   % ARmode��0���� instantaneous ģʽ������Ԫ���㣬1���� continuous ģʽ������Ԫ���㣬2���� fix and hold ģʽ���̶�����ģʽ
Mode=0;     % mode��1���� static ģʽ����̬��0���� kinematic ģʽ����̬
% Mode=0;    
% rtk.rb=[-2005009.8220 , 5411170.2092, 2707855.2068];%data��׼վλ��
rtk.rb=[-2005033.7825, 5411163.2342, 2707856.6419];%data2
% rtk.rb=[-1835468.425, 5599626.507, 2432544.873];
% rtk.rb=[-1991608.3225587 5418552.665654 2702743.468404];
% rtk.rb=[-2005039.0280, 5411161.1383, 2707856.5661];
for k=1:epochNUM   %��Ԫ��                                                                                                                                                                                                  
    
    if k==1
        pretime=[0,0];%��һ��ԪpretimeΪ0
    else  
        pretime=rtk.sol.time;%������ǵ�һ��Ԫ������ʱ����ǵ�ǰʱ�䣻
    end
    
     rtk.sol=pointpos(rtk.sol.rr,ObsODat_rover(k),EphDat,N);
% %      
%      xchange=rtk.sol(k).rr(1);
%      ychange=rtk.sol(k).rr(2);
%      zchange=rtk.sol(k).rr(3);   
%      figure(1);
% end
%       plot(xchange,'b.-');
%       xlabel('��Ԫ');ylabel('change(m)');
%       hold on;
%       plot(ychange,'r.-');
%       xlabel('��Ԫ');ylabel('change(m)');
%       hold on;
%       plot(zchange,'k.-');
%       xlabel('��Ԫ');ylabel('change(m)'); 
%       hold on;

%     ***********************
    if pretime(1)~=0
        rtk.tt=timediff(rtk.sol.time,pretime);%��ǰʱ����֮ǰʱ��֮��
    end

    rtk=relativepos(rtk,ObsODat_rover(k),ObsODat_base(k),EphDat,N,Mode,ARmode);

    
    for w=1:MAXSAT
        rtk.ssat(w).azel=[0,0];
    end
%     *************************
    sol(k)=rtk.sol;
    
    xchange(k)=Xk0-sol(k).rr(1);
    ychange(k)=Yk0-sol(k).rr(2);
    zchange(k)=Zk0-sol(k).rr(3);
    ratio(k)=rtk.sol.ratio;
    
end

   figure(1);

plot(xchange,'b.-');
xlabel('��Ԫ');ylabel('change(m)');
hold on

plot(ychange,'r.-');
xlabel('��Ԫ');ylabel('change(m)');
hold on

plot(zchange,'k.-');
xlabel('��Ԫ');ylabel('change(m)');
legend('xchange','ychange','zchange');
toc