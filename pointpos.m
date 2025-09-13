function sol=pointpos(solrr,ObsODat_rover_k,EphDat,N)
% *************************************************************************
% input:
% solrr                        上一历元的接收机位置，第一历元为0
% ObsODat_rover_k              本历元的接收机观测值
% EphDat                       导航数据
% N                            导航数据中的卫星个数（包括重复的）
% *************************************************************************
% output:
% sol                          单点定位的解
% *************************************************************************
%只能读取U-blox格式的单GPS数据

tr=ObsODat_rover_k.TimeOEp;                %历元接收数据时间（本历元，相当于obs[0].time）
Snum=ObsODat_rover_k.SatSum;               %该历元观测到的卫星数（N与Sat.Sum不同）
if Snum<4 
    return;                                %去除无解的历元
end
Code=ObsODat_rover_k.SatCode;              %该历元观测到的卫星号
R=ObsODat_rover_k.Obs_RangeC1;             %取出C1观测值,列向量（每种观测量都放在一列中）

time=epoch2time(tr(1),tr(2),tr(3),tr(4),tr(5),tr(6));%信号接收时间年月日时分秒
tJD=TimetoJD(tr(1),tr(2),tr(3),tr(4),tr(5),tr(6));%儒略日时间
sol.time=time;%时间更新为信号接收时间

%卫星位置
[rs,dts,vare]=satposs(sol.time,tJD,ObsODat_rover_k,Snum,EphDat,R,N);

%接收机位置
sol=estpos(time,solrr,rs,dts,Snum,R,Code,vare,EphDat);
