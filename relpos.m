%function relpos
% version history ***********************************************************
% 2017/05/12
% 为简化程序，使输入的移动站和基站数据开始与结束时间相同，历元数相同，直接双差
% 2017/05/24
% 增加了contineous模糊度解算模式，增加了static静态解算模式
% 2017/05/28
% 增加了fix and hold 模糊度解算模式
% 2017/06/05
% 修复了卡尔曼滤波中的bug，循环迭代中的卡尔曼滤波与holdamb函数中的卡尔曼滤波对矩阵H的求法不同
% **************************************************************************

close all;clear ;clc;
tic

 % if 0
global HeadODat_rover;%移动站的头文件
global ObsODat_rover; %移动站的观测数据
global HeadODat_base; %基站的头文件
global ObsODat_base;  %基站的观测数据
global EphDat;  %导航数据

global CLIGHT;%光速
global SOLQ_NONE;%无解 
global SOLQ_SINGLE;%单点解
global SOLQ_FLOAT;%浮点解
global SOLQ_FIX;%固定解
global MINFIX;%最小固定数
global MAXSAT;%最多可接收的卫星数 
global FREQ1;%L1频段
global FREQ2;%L2频段
global D2R;
global R2D;
global RE_WGS84;%地球长半径
global FE_WGS84;%扁率
global OMGE;%地球自转角速度大小

CLIGHT=299792458;%光速
SOLQ_NONE=0;
SOLQ_SINGLE=1;
SOLQ_FLOAT=2;
SOLQ_FIX=3;
MINFIX=10;
MAXSAT=32;
FREQ1=1.57542E+9;    %L1频段
FREQ2=1.22760E+9;    %L2频段
D2R=pi/180.0;%1度所对应的弧度
R2D=180.0/pi;%1弧度所对应的度
RE_WGS84=6378137.0;%长半径
FE_WGS84=(1.0/298.257223563);%扁率
OMGE=7.2921151467E-5;%地球自转角速度

EphDat=ReadGpsData;                                       %读N文件
[HeadODat_rover,ObsODat_rover]=ReadObsData(0);            %读移动站O文件
[HeadODat_base,ObsODat_base]=ReadObsData(1);              %读基站O文件

 % end
%  load ('F:\李-资料\matlab程序\rtklib_relpos\data\text.mat');
% load ('F:\李-资料\matlab程序\rtklib_relpos\data\text2.mat');
N = size(EphDat,2);                                       %得到EphDat(i)中i的个数（N文件里卫星总个数，在o文件中，因历元不同卫星数可能不同）
epochNUM=size(ObsODat_rover,2);                           %得到历元数，即观测数据的个数（每历元得到一次数据）

% Xk0=-2005024.6830;                                        %测站点的近似坐标
% Yk0=5411166.0497;
% Zk0=2707856.7714;

% Xk0=-2005009.6714;
% Yk0=5411172.1014;
% Zk0=2707856.3983; %data解锁位置

Xk0=-2005033.7825+2.1130;%enu
Yk0=5411163.2342+0.0233;
Zk0=2707856.6419+0.0242;
%该组数据为data2文件夹中的rinex2.10文件夹
% Xk0=-2005035.76699731;
% Yk0=5411162.50934514;
% Zk0=2707856.67218072;
%该组数据为data2文件夹中的a.nav和base，rover
% Xk0=-1991315.213867; 
% Yk0=5418624.208629; 
% Zk0=2702815.522514;



rtk=rtkinit;%初始化
ARmode=2;   % ARmode，0代表 instantaneous 模式即单历元解算，1代表 continuous 模式即多历元解算，2代表 fix and hold 模式即固定保持模式
Mode=0;     % mode，1代表 static 模式即静态，0代表 kinematic 模式即动态
% Mode=0;    
% rtk.rb=[-2005009.8220 , 5411170.2092, 2707855.2068];%data基准站位置
rtk.rb=[-2005033.7825, 5411163.2342, 2707856.6419];%data2
% rtk.rb=[-1835468.425, 5599626.507, 2432544.873];
% rtk.rb=[-1991608.3225587 5418552.665654 2702743.468404];
% rtk.rb=[-2005039.0280, 5411161.1383, 2707856.5661];
for k=1:epochNUM   %历元数                                                                                                                                                                                                  
    
    if k==1
        pretime=[0,0];%第一历元pretime为0
    else  
        pretime=rtk.sol.time;%如果不是第一历元，解算时间就是当前时间；
    end
    
     rtk.sol=pointpos(rtk.sol.rr,ObsODat_rover(k),EphDat,N);
% %      
%      xchange=rtk.sol(k).rr(1);
%      ychange=rtk.sol(k).rr(2);
%      zchange=rtk.sol(k).rr(3);   
%      figure(1);
% end
%       plot(xchange,'b.-');
%       xlabel('历元');ylabel('change(m)');
%       hold on;
%       plot(ychange,'r.-');
%       xlabel('历元');ylabel('change(m)');
%       hold on;
%       plot(zchange,'k.-');
%       xlabel('历元');ylabel('change(m)'); 
%       hold on;

%     ***********************
    if pretime(1)~=0
        rtk.tt=timediff(rtk.sol.time,pretime);%当前时间与之前时间之差
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
xlabel('历元');ylabel('change(m)');
hold on

plot(ychange,'r.-');
xlabel('历元');ylabel('change(m)');
hold on

plot(zchange,'k.-');
xlabel('历元');ylabel('change(m)');
legend('xchange','ychange','zchange');
toc