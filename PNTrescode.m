function [v,var,nv,H]=PNTrescode(X,t,rs,dts,vare,Snum,R,Code,EphDat)
% *************************************************************************
% X                            接收机位置
% t                            信号接收时间
% rs                           卫星位置
% dts                          卫星钟差
% vare                         卫星位置的方差
% Snum                         卫星数目
% R                            伪距观测值
% Code                         卫星号数组
% EphDat                       导航数据
% *************************************************************************
% output:
% v                            伪距残差
% var                          残差的方差
% nv                           残差数
% H                            系数矩阵（nv*4）
% *************************************************************************
%azel，数组按照方位角、仰角排列，大小为2*Snum

global CLIGHT;
global FREQ1;
global FREQ2;
global RE_WGS84;
global OMGE;
H=zeros(1,4);
rr(1)=X(1,1);rr(2)=X(2,1);rr(3)=X(3,1);
dtr=X(4,1);%接收机位置和钟差
nv=1;%残差数，初始为1

pos=ecef2pos(rr);%经纬高

for j=1:Snum
    azel(2*j-1)=0;
    azel(j*2)=0;
    
    %求星地距离r和星地矢量e
    for m=1:3
        e(m)=rs(j,m)-X(m,1);%星地矢量
    end
    r=norm(e);%星地距离
    for m=1:3
        e(m)=e(m)/r;%星地单位矢量
    end
    r=r+OMGE*((rs(j,1)*X(2,1))-(rs(j,2)*X(1,1)))/CLIGHT;%sagnac effect 根据效应距离校正
    if(r<=0)
        continue;%跳过后面语句回到for循环
    end
    
    %求仰角和方位角
    el=pi/2;
    az=0;
    if(pos(3)>-RE_WGS84)
        enu=ecef2enu(pos,e);%东北天坐标
        if(dot(enu,enu)<1E-12)
            az=0;
        else
            az=atan2(enu(1),enu(2));
        end
        if(az<0)
            az=az+2*pi;%方位角
        end
        el=asin(enu(3));%仰角
    end
    
    azel(2*j-1)=az;
    azel(j*2)=el;
    
    if(el<(15*pi/180.0))
        continue;
    end
    
    %码偏移纠正的伪距differential code bias，具体参看《GNSS卫星中DCB的使用方法_胡玉坤》
    lam(1)=CLIGHT/FREQ1;%L1波长
    lam(2)=CLIGHT/FREQ2;%L2波长
    gamma=(lam(2)*lam(2))/(lam(1)*lam(1));%(f1/f2)的平方
    P1=R(j);
    if(R(j)==0)
        continue;
    end
    
    %Number=1;
    
    for i=1:size(EphDat,2)%size函数1是返回EphDat行数，2是返回列数
        if(EphDat(i).PRN~=Code(j))
            continue;
        end
        tgd=CLIGHT*EphDat(i).tgd(1);%EphDat.tgd群波延时校正值
        break;%直接退出for循环
    end
    P1_P2=(1-gamma)*tgd;    %没有p1-p2的值，则用tgd代替
    PC=P1-P1_P2/(1.0-gamma);%纠正伪距
    vmeas=0.3*0.3;          %伪距方差
    
    %电离层纠正

        ion_default=[0.1118E-07 -0.7451E-08 -0.5961E-07 0.1192E-06 0.1167E+06 -0.2294E+06 -0.1311E+06 0.1049E+07];
        if(pos(3)<(-1E+3)||azel(j*2)<=0)
            dion=0;
        else
            ion=ion_default;
            psi=0.0137/(azel(2*j)/pi+0.11)-0.022;
            phi=pos(1)/pi+psi*cos(azel(2*j-1));
            if(phi>0.416)
                phi=0.416;
            elseif(phi<-0.416)
                phi=-0.416;
            end
            LAM=pos(2)/pi+psi*sin(azel(2*j-1))/cos(phi*pi);
            phi=phi+0.064*cos((LAM-1.617)*pi);
            week=time2gpst(t);%用到输入量，信号接收时间
            tt=43200.0*LAM+week(2);
            tt=tt-floor(tt/86400.0)*86400.0;%0<=tt<86400.0
            f=1.0+16.0*(0.53-azel(j*2)/pi)^3;
            amp=ion(1)+phi*(ion(2)+phi*(ion(3)+phi*ion(4)));
            per=ion(5)+phi*(ion(6)+phi*(ion(7)+phi*ion(8)));
            if(amp<0)
                amp=0;
            end
            if(per<72000.0)
                per=72000.0;
            end
            x=2.0*pi*(tt-50400)/per;
            if(abs(x)<1.57)
                temp=5E-9+amp*(1.0+x*x*(-0.5+x*x/24.0));
            else
                temp=5E-9;
            end
            dion=CLIGHT*f*temp;%电离层修正
        end
        vion=(dion*0.5)*(dion*0.5);%电离层方差
    
    %GPS-L1 -> L1/B1
    lam_L1=lam(1);
    lam_carr=CLIGHT/FREQ1;
    if(lam_L1>0.0)
        dion=dion*(lam_L1/lam_carr)*(lam_L1/lam_carr);
    end
    
    %对流层修正

        if(pos(3)<-100||pos(3)>1E+4||azel(j*2)<=0)
            dtrp=0;
        else
            if(pos(3)<0)
                hgt=0;
            else
                hgt=pos(3);
            end
            pres=1013.25*(1.0-2.2557E-5*hgt)^5.2568;
            temp=15.0-6.5E-3*hgt+273.16;
            Edtrp=6.108*0.7*exp((17.15*temp-4684.0)/(temp-38.45));
            z=pi/2.0-azel(j*2);
            trph=0.0022768*pres/(1.0-0.00266*cos(2.0*pos(1))-0.00028*hgt/1E+3)/cos(z);
            trpw=0.002277*(1255.0/temp+0.05)*Edtrp/cos(z);
            dtrp=trph+trpw;%对流层修正
        end
        vtrp=(0.3/(sin(azel(j*2))+0.1))*(0.3/(sin(azel(j*2))+0.1));
    
    %伪距残差   残差都是用来算误差的
    v(nv)=PC-(r+dtr-CLIGHT*dts(j)+dion+dtrp);%dts输入量
    
%     if j==1
%         H=[-e(1) -e(2) -e(3) 1];
%     else
%         H=[H; -e(1) -e(2) -e(3) 1];        
%     end
%%%%%%%%%%改为下面这个2018-1-30
    if H~=0
        H=[H; -e(1) -e(2) -e(3) 1];
    else
        H=[-e(1) -e(2) -e(3) 1];
    end
    
    var(nv)=(100*100)*((0.003*0.003)+(0.003*0.003)/sin(azel(j*2)))+vare(j)+vmeas+vion+vtrp;
    nv=nv+1;
    %Number=Number+1;
end
nv=nv-1;
