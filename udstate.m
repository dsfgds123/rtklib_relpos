function rtk=udstate(mode,ARmode,rtk,ObsODat_rover_k,ObsODat_base_k,sat,ir,ib,ns,EphDat)

global FREQ1;
global CLIGHT;
global MAXSAT;

VAR_POS=900.0;      %(卫星位置方差)

rtkX=rtk.x;%接收机状态
rtkP=rtk.P;%方差协方差阵

tt=abs(rtk.tt);%当前时间与之前时间差

%%
%(1)更新位置、速度和加速度的状态和方差协方差矩阵
Norm=sqrt(rtk.x(1)*rtk.x(1)+rtk.x(2)*rtk.x(2)+rtk.x(3)*rtk.x(3));%接收机状态？？？这是判断是不是首历元

%initialize position for first epoch
if Norm<=0                                         %这都开了平方根了，怎么会小于0？？
    for i=1:3
        [rtkX,rtkP]=initx(rtk.nx,rtkX,rtkP,rtk.sol.rr(i),VAR_POS,i);%初始化
    end
else
    if mode==1                    %静态直接用上一历元的浮点解
        rtkX=rtk.x;
        rtkP=rtk.P;
    else                          %如果不是静态，且不是首历元，也需要初始化
        for i=1:3
            [rtkX,rtkP]=initx(rtk.nx,rtkX,rtkP,rtk.sol.rr(i),VAR_POS,i);
        end
    end
end
%%

%%
%(2)相位偏差更新
for i=1:ns
    %detect cycle slip by LLI  周跳检测
    rtk.ssat(sat(i)).slip=bitand(rtk.ssat(sat(i)).slip, 252);
    rtk.ssat(sat(i)).slip=detslp_ll(rtk,ObsODat_rover_k,ir(i),1);
    rtk.ssat(sat(i)).slip=detslp_ll(rtk,ObsODat_base_k,ib(i),2); 
end

%reset phase-bias if instantaneous AR or expire obs outage counter 
for i=1:MAXSAT
    %opt参数，5是obs outage counter of phase
    rtk.ssat(i).outc=rtk.ssat(i).outc+1;%outc卫星断供计数
    if (rtk.ssat(i).outc)>5
        reset=1;
    else
        reset=0;
    end
    
 %   if (ARmode==0&rtkX(3+i)~=0)
    if (ARmode==0&rtkX(3+i)~=0)
        [rtkX,rtkP]=initx(rtk.nx,rtkX,rtkP,0,0,3+i);
    elseif ((reset==1)&rtkX(3+i)~=0)
        [rtkX,rtkP]=initx(rtk.nx,rtkX,rtkP,0,0,3+i);
    end
    if ((ARmode~=0)&(reset==1))
        %min lock count to fix ambiguity,=0
        rtk.ssat(i).lock=-0;%载波相位锁定计数
    end
end

%reset phase-bias if detecting cycle slip 
for i=1:ns
    j=3+sat(i);
    %rtk->opt.prn[0],=1e-4
    rtkP(j,j)=rtkP(j,j)+1e-4*1e-4*tt;
    slip=rtk.ssat(sat(i)).slip;
    if ARmode==0|~(slip&1)%或语句，要么ARmode=0,要么slip=0，满足其一即可往下进行
        continue;
    end
%     if ARmode==0|~(slip&1)
%         continue;
%     end
    rtkX(j)=0;
    rtk.ssat(sat(i)).lock=-0;%？
end

bias=zeros(ns,1);

%estimate approximate phase-bias by phase - code 
offset=0;%offset是本历元和上一历元相位偏差的差值
j=0;
for i=1:ns
    cp=sgobs(ObsODat_rover_k,ObsODat_base_k,ir(i),ib(i),0);%载波相位差值
    pr=sgobs(ObsODat_rover_k,ObsODat_base_k,ir(i),ib(i),1);%伪距差值
    lami=CLIGHT/FREQ1;%波长
    if cp==0|pr==0|lami<=0
        continue;
    end
    bias(i)=cp-pr/lami;
    
    if rtkX(3+sat(i))~=0
        offset=offset+bias(i)-rtkX(3+sat(i));
        j=j+1;
    end
end

%correct phase-bias offset to enssure phase-code coherency
if j>0
    for i=1:MAXSAT
        if rtkX(3+i)~=0
            rtkX(3+i)=rtkX(3+i)+offset/j;%更新结果为bias+(offset的均值)
        end
    end
end

%set initial states of phase-bias
for i=1:ns
    if bias(i)==0|rtkX(3+sat(i))~=0
        continue;
    end
    %SQR(rtk->opt.std[0]),=900
    [rtkX,rtkP]=initx(rtk.nx,rtkX,rtkP,bias(i),900,3+sat(i));
end

rtk.x=rtkX;
rtk.P=rtkP;
%%
%%
function cp_pr=sgobs(ObsODat_rover_k_,ObsODat_base_k_,ir_,ib_,zeroORone)

if zeroORone==0
    pi=ObsODat_rover_k_.Obs_FreL1(ir_);
    pj=ObsODat_base_k_.Obs_FreL1(ib_);
elseif zeroORone==1
    pi=ObsODat_rover_k_.Obs_RangeC1(ir_);
    pj=ObsODat_base_k_.Obs_RangeC1(ib_);
end

if pi==0|pj==0
    cp_pr=0;
else
    cp_pr=pi-pj;
end
%%

