function [nv,rtk,v,H,R,vflg]=ddres(rtk,dt,x,sat,y_r,e_r,azel_r,y_b,ir,ib,ns,Hmode)
% *************************************************************************
% input:
% rtk                          解算参数结构体
% dt                           基站移动站接收数据数据时间差
% x                            移动站接收机位置
% sat                          公式卫星的卫星号
% y_r                          移动站接收机非差分残差（包括伪距和载波相位）
% e_r                          移动站接收机与卫星之间的向量
% azel_r                       移动站接收机与卫星之间的高度角和仰角
% y_b                          基站接收机非差分残差
% ir                           共视卫星的移动站接收机卫星号排序
% ib                           共视卫星的基站接收机卫星号排序
% ns                           共视卫星个数
% Hmode                        标记求或不求系数矩阵H
% *************************************************************************
% output:
% nv                           双差残差个数
% rtk                          解算参数结构体
% v                            双插残差向量
% H                            系数矩阵H
% R                            双差测量误差的方差协方差
% vflg                         标志位
% *************************************************************************

global CLIGHT;
global MAXSAT;
global FREQ1;

H=zeros(ns*2,rtk.nx);%%共视卫星数目*2，状态数量，生成一个（18，35）的零矩阵

didxi=0;
didxj=0;
nv=1;
b=1;     %%？？
nb=[0,0];%%？？

[bl,dr]=baseline(x,rtk.rb);%x移动站位置和rb基站位置  bl基线距离   dr基线矢量

posu=ecef2pos(x);
posr=ecef2pos(rtk.rb);%转为经纬高

for i=1:MAXSAT
    rtk.ssat(i).resp=0;          %rtk.ssat为卫星状态结构体，rtk.ssat.resp伪距双差残差，rtk.ssat.resc 载波相位双差残差
    rtk.ssat(i).resc=0;
end

for f=1:2
    %search reference satellite with highest elevation
    %将海拔最高的卫星当做参考卫星即仰角最高的共视卫星
    i=-1;
    for j=1:ns
        if i<0|azel_r(2*ir(j))>=azel_r(2*ir(i))
            i=j;
        end
    end    %确定i的值即确定仰角最高的共视卫星（这里所得到的仰角最高卫星i=3）
    if i<0
        continue;
    end
    
    %make double difference
    for j=1:ns
        if i==j
            continue;
        end

        lami=CLIGHT/FREQ1;  %%波长
        lamj=CLIGHT/FREQ1;
        
        %double-differenced residual%   f=1 得到载波相位双差残差，f=2得到伪距双差残差
        v(nv)=(y_r(f+(ir(i)-1)*2)-y_b(f+(ib(i)-1)*2))...   %这里是先做站间差再做星间差
            -(y_r(f+(ir(j)-1)*2)-y_b(f+(ib(j)-1)*2)); 
        
        %partial derivatives by rover position 移动站位置偏导
        if Hmode==1
            for k=1:3
                H(nv,k)=-e_r(ir(i),k)+e_r(ir(j),k);%nv为残差数关于H矩阵看文档有详细的计算
            end
        end
        
        %double-differenced phase-bias term 双差相位残差项，对模糊度的求解
        if f==1
            v(nv)=v(nv)-lami*x(3+sat(i))+lami*x(3+sat(j));%  测量值 减 预测值   后面两项是个差值，我觉得是给双差残差修正
            if Hmode==1
                H(nv,3+sat(i))= lami;
                H(nv,3+sat(j))=-lamj;
            end
        end
        
        if f==1
            rtk.ssat(sat(j)).resc=v(nv);% resc载波相位残差
        else
            rtk.ssat(sat(j)).resp=v(nv);% resp伪距残差 更新RTK结构体
        end
        
        %test innovation
        % reject threshold of innovation, =30
        if abs(v(nv))>30
            if f==1
                rtk.ssat(sat(i)).rejc=rtk.ssat(sat(i)).rejc+1;%排除的数据的个数
                rtk.ssat(sat(j)).rejc=rtk.ssat(sat(j)).rejc+1;
            end
%             continue;
        end
        
        %single-differenced measurement error variances 单差测量误差方差
        Ri(nv)=vareer(azel_r(ir(i)*2),bl,dt,f);%dt基站和移动站接收数据的时间差  仰角
        Rj(nv)=vareer(azel_r(ir(j)*2),bl,dt,f);
        
        %set valid data flags 设置有效数据标志
        if f==1
            rtk.ssat(sat(i)).vsat=1;
            rtk.ssat(sat(j)).vsat=1;
        end
        
        if f==1
            w=0;
        else
            w=1;
        end
        vflg(nv)=bitor(bitor(bitor(bitshift(sat(i),16),bitshift(sat(j),8)),bitshift(w,4)),0);
        nv=nv+1;
        nb(b)=nb(b)+1;
    end
    b=b+1;
end
b=b-1;
nv=nv-1;

%double-differenced measurement error covariance
R=ddcov(nb,b,Ri,Rj,nv);


function [bl_,dr_]=baseline(rr,rb)
for i=1:3
    dr_(i)=rr(i)-rb(i);%向量
end
bl_=norm(dr_);%返回最大奇异值基站和移动站的直线距离（基线）


function Ri_nv=vareer(el,bl,dt,f)  %单差误差协方差函数
CLIGHT=299792458.0;
c=0*bl/1e4;%c=0
d=CLIGHT*5e-12*dt;
fact=1;

if f==2
    fact=100;
end
if fact<=0
    fact=100;
end

%opt->err[1]=opt->err[2]=0.003
a=fact*0.003;
b=fact*0.003;

Ri_nv=2*(a*a+b*b/sin(el)/sin(el)+c*c)+d*d;%单差误差协方差阵

function R_=ddcov(nb_,b_,Ri_,Rj_,nv_)  %双差误差协方差函数
k=0;
R_=zeros(nv_);

for b=1:b_
    for i=1:nb_(b_)
        for j=1:nb_(b_)
            if i==j
                R_(k+i,k+j)=Ri_(k+i)+Rj_(k+i);
            else
                R_(k+i,k+j)=Ri_(k+i);
            end
        end
    end
    k=k+nb_(b);
end
