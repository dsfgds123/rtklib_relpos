function rtk=udstate(mode,ARmode,rtk,ObsODat_rover_k,ObsODat_base_k,sat,ir,ib,ns,EphDat)

global FREQ1;
global CLIGHT;
global MAXSAT;

VAR_POS=900.0;      %(����λ�÷���)

rtkX=rtk.x;%���ջ�״̬
rtkP=rtk.P;%����Э������

tt=abs(rtk.tt);%��ǰʱ����֮ǰʱ���

%%
%(1)����λ�á��ٶȺͼ��ٶȵ�״̬�ͷ���Э�������
Norm=sqrt(rtk.x(1)*rtk.x(1)+rtk.x(2)*rtk.x(2)+rtk.x(3)*rtk.x(3));%���ջ�״̬�����������ж��ǲ�������Ԫ

%initialize position for first epoch
if Norm<=0                                         %�ⶼ����ƽ�����ˣ���ô��С��0����
    for i=1:3
        [rtkX,rtkP]=initx(rtk.nx,rtkX,rtkP,rtk.sol.rr(i),VAR_POS,i);%��ʼ��
    end
else
    if mode==1                    %��ֱ̬������һ��Ԫ�ĸ����
        rtkX=rtk.x;
        rtkP=rtk.P;
    else                          %������Ǿ�̬���Ҳ�������Ԫ��Ҳ��Ҫ��ʼ��
        for i=1:3
            [rtkX,rtkP]=initx(rtk.nx,rtkX,rtkP,rtk.sol.rr(i),VAR_POS,i);
        end
    end
end
%%

%%
%(2)��λƫ�����
for i=1:ns
    %detect cycle slip by LLI  �������
    rtk.ssat(sat(i)).slip=bitand(rtk.ssat(sat(i)).slip, 252);
    rtk.ssat(sat(i)).slip=detslp_ll(rtk,ObsODat_rover_k,ir(i),1);
    rtk.ssat(sat(i)).slip=detslp_ll(rtk,ObsODat_base_k,ib(i),2); 
end

%reset phase-bias if instantaneous AR or expire obs outage counter 
for i=1:MAXSAT
    %opt������5��obs outage counter of phase
    rtk.ssat(i).outc=rtk.ssat(i).outc+1;%outc���ǶϹ�����
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
        rtk.ssat(i).lock=-0;%�ز���λ��������
    end
end

%reset phase-bias if detecting cycle slip 
for i=1:ns
    j=3+sat(i);
    %rtk->opt.prn[0],=1e-4
    rtkP(j,j)=rtkP(j,j)+1e-4*1e-4*tt;
    slip=rtk.ssat(sat(i)).slip;
    if ARmode==0|~(slip&1)%����䣬ҪôARmode=0,Ҫôslip=0��������һ�������½���
        continue;
    end
%     if ARmode==0|~(slip&1)
%         continue;
%     end
    rtkX(j)=0;
    rtk.ssat(sat(i)).lock=-0;%��
end

bias=zeros(ns,1);

%estimate approximate phase-bias by phase - code 
offset=0;%offset�Ǳ���Ԫ����һ��Ԫ��λƫ��Ĳ�ֵ
j=0;
for i=1:ns
    cp=sgobs(ObsODat_rover_k,ObsODat_base_k,ir(i),ib(i),0);%�ز���λ��ֵ
    pr=sgobs(ObsODat_rover_k,ObsODat_base_k,ir(i),ib(i),1);%α���ֵ
    lami=CLIGHT/FREQ1;%����
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
            rtkX(3+i)=rtkX(3+i)+offset/j;%���½��Ϊbias+(offset�ľ�ֵ)
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

