function pntpos
close all;clear all;clc;
tic
  global HeadODat;
  global ObsODat; 
  global EphDat;

  sol=struct;

EphDat=ReadGpsData;%先读N文件，再读O文件
[HeadODat,ObsODat]=ReadObsData(0);

%多个历元加权平均求测站点坐标

N = size(EphDat,2);%得到EphDat(i)中i的个数（N文件里卫星总个数，在o文件中，因历元不同卫星数可能不同）
CLIGHT=299792458;%光速
epochNUM=size(ObsODat,2);%得到历元数，即观测数据的个数（每历元得到一次数据）

% Xk0=HeadODat.ApproXYZ(1,1);%测站点的近似坐标
% Yk0=HeadODat.ApproXYZ(1,2);
% Zk0=HeadODat.ApproXYZ(1,3);

Xk0=-2005009.6714;
Yk0=5411172.1014;
Zk0=2707856.3983;

for k=1:epochNUM

    tr=ObsODat(k).TimeOEp;         %历元接收数据时间（本历元，相当于obs[0].time）
    Snum=ObsODat(k).SatSum;      %该历元观测到的卫星数（N与Sat.Sum不同）
    if Snum<4 
        break;                     %去除无解的历元
    end     
    Code=ObsODat(k).SatCode;       %该历元观测到的卫星号
    R=ObsODat(k).Obs_RangeC1 ;     %取出C1观测值,列向量（每种观测量都放在一列中）

    time=epoch2time(tr(1),tr(2),tr(3),tr(4),tr(5),tr(6));%信号接收时间
    tJD=TimetoJD(tr(1),tr(2),tr(3),tr(4),tr(5),tr(6));

    %卫星位置
    [rs,dts,vare]=satposs(time,tJD,ObsODat(k),Snum,EphDat,R,N);

    %接收机位置
    if k==1
        sol0rr=[0 0 0];
        sol=estpos(time,sol0rr,rs,dts,Snum,R,Code,vare,EphDat);
    else
        sol=estpos(time,sol(k-1).rr,rs,dts,Snum,R,Code,vare,EphDat);
    end

    sol(k)=sol;

    xchange(k)=Xk0-sol(k).rr(1);
    ychange(k)=Yk0-sol(k).rr(2);
    zchange(k)=Zk0-sol(k).rr(3);

end

figure(1);
subplot(3,1,1);plot(xchange,'b.-');
subplot(3,1,2);plot(ychange,'r.-');
subplot(3,1,3);plot(zchange,'k.-');
toc


