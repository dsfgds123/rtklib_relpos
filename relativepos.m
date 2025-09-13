function rtk=relativepos(rtk,ObsODat_rover_k,ObsODat_base_k,EphDat,N,mode,ARmode)

global MINFIX;
global SOLQ_NONE;%无解
global SOLQ_FLOAT;%浮点解
global SOLQ_FIX;%整数解
global MAXSAT;

stat=SOLQ_FLOAT;

Snum_r=ObsODat_rover_k.SatSum; Snum_b=ObsODat_base_k.SatSum;                 %当前历元基站移动站的卫星数

R_r=ObsODat_rover_k.Obs_RangeC1; R_b=ObsODat_base_k.Obs_RangeC1;             %基站移动站各卫星的伪距测量值

tr_r=ObsODat_rover_k.TimeOEp;                                                  %基站移动站本历元时间
time_r=epoch2time(tr_r(1),tr_r(2),tr_r(3),tr_r(4),tr_r(5),tr_r(6));
tr_b=ObsODat_base_k.TimeOEp;
time_b=epoch2time(tr_b(1),tr_b(2),tr_b(3),tr_b(4),tr_b(5),tr_b(6));

dt=timediff(time_r,time_b);                                                    %基站移动站接收数据数据时间差

tJD_r=TimetoJD(tr_r(1),tr_r(2),tr_r(3),tr_r(4),tr_r(5),tr_r(6));               %儒略日时间

[RS_r,DTS_r,VAR_r]=satposs(time_r,tJD_r,ObsODat_rover_k,Snum_r,EphDat,R_r,N); %移动站卫星位置  
[RS_b,DTS_b,VAR_b]=satposs(time_r,tJD_r,ObsODat_base_k,Snum_b,EphDat,R_b,N);  %基站卫星位置-
 %基站非差分残差
[y_b,e_b,azel_b]=zdres(time_b,ObsODat_base_k,Snum_b,RS_b,DTS_b,rtk.rb); %第一次的非差分残差是为了仰角方位角去寻找共视卫星                 

%time-interpolation of residuals (for post-processing)
%dt=intpres(time_r,ObsODat_base_k,Snum_b,EphDat,rtk,y_b);

[ns,sat,ir,ib]=selsat(ObsODat_rover_k,Snum_r,ObsODat_base_k,Snum_b,azel_b);%查找共视卫星

rtk=udstate(mode,ARmode,rtk,ObsODat_rover_k,ObsODat_base_k,sat,ir,ib,ns,EphDat);%更新RTK结构体（当前历元初始化）

% if abs(2005009.6688215812+rtk.x(1))>0.1 %做一个判断：应该是和真实坐标做一个比较
%     aaa=1;
%     bbb=2;   %？？？？ 
% end

xp=rtk.x;%接收机状态
% 1是迭代次数，可根据需要调整迭代次数
% for i=1:1                                    %这个循环是啥意思？？
    %移动站非差分残差
    [y_r,e_r,azel_r]=zdres(time_r,ObsODat_rover_k,Snum_r,RS_r,DTS_r,xp);%输出载波相位伪距残差，星地单位矢量以及仰角方位角
    
    for j=1:Snum_r
        for k=1:MAXSAT
            if ObsODat_rover_k.SatCode(j)==k
                rtk.ssat(k).azel=[azel_r(2*j-1),azel_r(2*j)];%？？？方位角  把解算出来的方位角仰角数据存入RTK里
                break;
            end
        end
    end

    %最后的参数1代表ddres中给H赋值，0代表不赋值
    [nv,rtk,v,H,R,vflg]=ddres(rtk,dt,xp,sat,y_r,e_r,azel_r,y_b,ir,ib,ns,1);%%双差计算
    
    Pp=rtk.P;%%方差协方差阵
    
    [xp_,Pp_]=EKF_filter(xp,Pp,H,v,R,rtk.nx,nv,0);%卡尔曼滤波，滤除后的位置与接收机方差
% end

if(stat~=SOLQ_NONE) %如果浮点解不是无解                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                  
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                    
    [y_r,e_r,azel_r]=zdres(time_r,ObsODat_rover_k,Snum_r,RS_r,DTS_r,xp_);%卡尔曼滤波过后再进行非差解算
    
    %此处H的值，与EKF_filter得到值相同，没有变化
    [nv,rtk,v,H,R,vflg]=ddres(rtk,dt,xp_,sat,y_r,e_r,azel_r,y_b,ir,ib,ns,0);%滤波之后的残差和双差
    
    if(valpos(rtk,v,R,vflg,nv,4))%有效性检验
        %更新状态和协方差矩阵
        rtk.x=xp_;
        rtk.P=Pp_;
        %更新模糊度控制结构体
        rtk.sol.ns=0;%有效卫星数
        for i=1:ns             %%共视卫星数
            if ~rtk.ssat(sat(i)).vsat%%数据是否有效标志
                continue;
            end
            rtk.ssat(sat(i)).lock=rtk.ssat(sat(i)).lock+1;%相位锁定计数
            rtk.ssat(sat(i)).outc=0;%卫星信号断供技术
            rtk.sol.ns=rtk.sol.ns+1;
        end
        %缺少有效卫星
        if rtk.sol.ns<4
            stat=SOLQ_NONE;
        end
    else
        stat=SOLQ_NONE;
    end
end

if stat~=SOLQ_NONE
    % [rtk,xa,AmbNum]=resamb_LAMBDA(rtk);
    [rtk,xa,AmbNum]=resamb_AGGWO(rtk);%求解整周模糊度及对应接收机状态
    if AmbNum>1
        [y_r,e_r,azel_r]=zdres(time_r,ObsODat_rover_k,Snum_r,RS_r,DTS_r,xa);%更新后的xa计算残差
        
        %这里也不求矩阵H，可以在ddres函数的传递参数中添加一个参数H_YorN来标记求或不求H
        %此处也不求用带入Pp，可以添加一个参数Pp_YorN
        [nv,rtk,v,H,R,vflg]=ddres(rtk,dt,xa,sat,y_r,e_r,azel_r,y_b,ir,ib,ns,0);%更新后的数值计算双差
        
        if(valpos(rtk,v,R,vflg,nv,4))%有效性检验
            %保持整周模糊度
            rtk.nfix=rtk.nfix+1;%固定解个数加1
            if(rtk.nfix>=MINFIX&ARmode==2)
                [SsatFix,rtk.x,rtk.P]=holdamb(rtk,xa);
                for i=1:MAXSAT
                    rtk.ssat(i).fix=SsatFix(i);%卫星状态结构体固定解
                end
            end
            stat=SOLQ_FIX;
        end
    end
end
%保存解
if stat==SOLQ_FIX
    rtk.sol.rr=[rtk.xa(1), rtk.xa(2), rtk.xa(3)];
    rtk.sol.qr=[rtk.Pa(1,1), rtk.Pa(2,2), rtk.Pa(3,3)];
else
    rtk.sol.rr=[rtk.x(1), rtk.x(2), rtk.x(3)];
    rtk.sol.qr=[rtk.P(1,1), rtk.P(2,2), rtk.P(3,3)];
    rtk.nfix=0;
end

for i=1:MAXSAT
    if rtk.ssat(i).fix==2&stat~=SOLQ_FIX       %如果卫星状态固定数为2，且浮点解不固定
        rtk.ssat(i).fix=1;                     %浮点解为1；
    end
end
