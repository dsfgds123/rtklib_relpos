function [RS,DTS,VAR]=satposs(time,tJD,ObsODatn,Snum,EphDat,R,N)%读卫星位置的函数
% input:
% time                         本历元的信号接收时间
% tJD                          信号接收时间的儒略日
% ObsODatn                     第n个历元的观测数据
% Snum                         本历元观测到的卫星个数
% EphDat                       导航数据
% R                            L1的伪距观测值
% N                            导航数据的卫星个数
% *************************************************************************
% output:
% RS                           卫星位置
% DTS                          卫星钟差
% VAR                          卫星位置的方差
% *************************************************************************

global CLIGHT;%光速

Code=ObsODatn.SatCode;       %该历元观测到的卫星号

for j=1:Snum
    
    for i=1:N               %每个历元所观测到的卫星还必须是导航数据文件的卫星
        timei= timeadd(time,-R(j)/CLIGHT);%transmission time by satellite clock ，得到计算的信号发射时间
        if(EphDat(i).PRN==Code(j)&&abs(tJD-EphDat(i).toc)<0.0417*2)%读入时间相近的星历文件
            %satellite clock bias by broadcast ephemeris  由广播星历表引起的卫星时钟偏差
            a0= EphDat(i).a0;
            a1= EphDat(i).a1;
            a2= EphDat(i).a2;
            %toe=EphDat(i).toe;
            ta=((timei(1)-EphDat(i).tocG(1))+(timei(2)-EphDat(i).tocG(2)));%toc钟差参数的参考时间，t-toc  toe 星历参考时间
            for w=1:2
                ta=ta-(a0+a1*ta+a2*ta*ta);  %多次修正卫星钟差 (a0+a1*ta+a2*ta*ta)为卫星钟差
            end
%             tb=ta-(a0+a1*ta+a2*ta*ta);
%             tc=tb-(a0+a1*tb+a2*tb*tb);
            dt=a0+ta*a1+ta*ta*a2;%卫星钟差
            
            timei=timeadd(timei,-dt);
            
            %satellite position and clock at transmission time 卫星位置和时钟在传输时间
            [rs(1),rs(2),rs(3),dts,var]= satpos(EphDat(i),timei);
            if j==1
                RS=[rs(1),rs(2),rs(3)];
            else
                RS=[RS;rs(1),rs(2),rs(3)];
            end
            DTS(j)=dts;
            VAR(j)=var;
            %ts=ts2+tm;%以下为求钟漂
            %[Xs2,Ys2,Zs2,dts2]= CalPos1(EphDat(i),ts2);
            %Xs=(Xs2-Xs)/tm; %Ys=(Ys2-Ys)/tm; %Zs=(Zs2-Zs)/tm;
            %dts=(dts2-dts)/tm;
            %break;
        end
    end
end
