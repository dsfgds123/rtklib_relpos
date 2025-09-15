function pntpos
close all;clear all;clc;
tic
  global HeadODat;
  global ObsODat; 
  global EphDat;

  sol=struct;

EphDat=ReadGpsData;%�ȶ�N�ļ����ٶ�O�ļ�
[HeadODat,ObsODat]=ReadObsData(0);

%�����Ԫ��Ȩƽ�����վ������

N = size(EphDat,2);%�õ�EphDat(i)��i�ĸ�����N�ļ��������ܸ�������o�ļ��У�����Ԫ��ͬ���������ܲ�ͬ��
CLIGHT=299792458;%����
epochNUM=size(ObsODat,2);%�õ���Ԫ�������۲����ݵĸ�����ÿ��Ԫ�õ�һ�����ݣ�

% Xk0=HeadODat.ApproXYZ(1,1);%��վ��Ľ�������
% Yk0=HeadODat.ApproXYZ(1,2);
% Zk0=HeadODat.ApproXYZ(1,3);

Xk0=-2005009.6714;
Yk0=5411172.1014;
Zk0=2707856.3983;

for k=1:epochNUM

    tr=ObsODat(k).TimeOEp;         %��Ԫ��������ʱ�䣨����Ԫ���൱��obs[0].time��
    Snum=ObsODat(k).SatSum;      %����Ԫ�۲⵽����������N��Sat.Sum��ͬ��
    if Snum<4 
        break;                     %ȥ���޽����Ԫ
    end     
    Code=ObsODat(k).SatCode;       %����Ԫ�۲⵽�����Ǻ�
    R=ObsODat(k).Obs_RangeC1 ;     %ȡ��C1�۲�ֵ,��������ÿ�ֹ۲���������һ���У�

    time=epoch2time(tr(1),tr(2),tr(3),tr(4),tr(5),tr(6));%�źŽ���ʱ��
    tJD=TimetoJD(tr(1),tr(2),tr(3),tr(4),tr(5),tr(6));

    %����λ��
    [rs,dts,vare]=satposs(time,tJD,ObsODat(k),Snum,EphDat,R,N);

    %���ջ�λ��
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


