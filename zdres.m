function [y,e,azel]=zdres(time,ObsODat_k,n,RS,DTS,rr)
% *************************************************************************
% input:
% time                 时间信息
% base                标志输入的是基站还是移动站数据（0，移动站；1，基站）
% HeadODatk           本观测历元的观测数据
% n                   观测数据的卫星数
% RS                  卫星位置
% DTS                 卫星钟差
% rr                  接收机位置
% *************************************************************************
% output:
% y                   伪距和载波相位非差分残差
% e                   星地矢量
% azel                仰角方位角
% *************************************************************************

global CLIGHT;
global FREQ1;
global D2R;
global R2D;
global RE_WGS84;
global OMGE;

coef=[1.2769934E-3, 1.2683230E-3, 1.2465397E-3, 1.2196049E-3, 1.2045996E-3;
      2.9153695E-3, 2.9152299E-3, 2.9288445E-3, 2.9022565E-3, 2.9024912E-3;
      62.610505E-3, 62.837393E-3, 63.721774E-3, 63.824265E-3, 64.258455E-3;
      0.0000000E-0, 1.2709626E-5, 2.6523662E-5, 3.4000452E-5, 4.1202191E-5;
      0.0000000E-0, 2.1414979E-5, 3.0160779E-5, 7.2562722E-5, 11.723375E-5;
      0.0000000E-0, 9.0128400E-5, 4.3497037E-5, 84.795348E-5, 170.37206E-5;
      5.8021897E-4, 5.6794847E-4, 5.8118019E-4, 5.9727542E-4, 6.1641693E-4;
      1.4275268E-3, 1.5138625E-3, 1.4572752E-3, 1.5007428E-3, 1.7599082E-3;
      4.3472961E-2, 4.6729510E-2, 4.3908931E-2, 4.4626982E-2, 5.4736038E-2];
aht=[2.53E-5, 5.49E-3, 1.14E-3];

zazel=[0,90*D2R];
rr_=[rr(1), rr(2), rr(3)];

%earth tide correction
%tidedisp()

pos=ecef2pos(rr_);

for i=1:n
    
    %求星地距离和星地矢量
    for m=1:3
        e(i,m)=RS(i,m)-rr_(m);%%星站矢量（x，y，z）
    end
    r=norm(e(i,:));%%e(i;:)表示e中的i行全部数据。norm星地距离
    for m=1:3
        e(i,m)=e(i,m)/r;%%矢量除以距离，为单位矢量
    end
    r=r+OMGE*((RS(i,1)*rr_(2))-(RS(i,2)*rr_(1)))/CLIGHT;%%地球自转修正
    if(r<=0)
        continue;
    end
    
    %求仰角方位角 
    el=pi/2;
    az=0;
    if(pos(3)>-RE_WGS84)
        enu=ecef2enu(pos,e(i,:));
        if((enu(1)*enu(1)+enu(2)*enu(2))<1E-12)
            az=0;
        else
            az=atan2(enu(1),enu(2)); 
        end
        if(az<0)
            az=az+2*pi;
        end
        el=asin(enu(3));
    end
    azel(2*i-1)=az;
    azel(i*2)=el;
    
    if(el<(15*pi/180.0))
        continue;
    end
    
    %satellite clock-bias修正星地距离
    r=r-CLIGHT*DTS(i);%修正卫星钟差
    
    %对流层延迟模型修正星地距离(hydrostatic)
    if(pos(3)<-100||pos(3)>1E+4||zazel(2)<=0)
%         dtrp=0;???
        zhd = 0;
    else
        if(pos(3)<0)
            hgt=0;
        else
            hgt=pos(3);
        end
        pres=1013.25*(1.0-2.2557E-5*hgt)^5.2568;
        temp=15.0-6.5E-3*hgt+273.16;
        Edtrp=6.108*0*exp((17.15*temp-4684.0)/(temp-38.45));
        z=pi/2.0-zazel(2);
        trph=0.0022768*pres/(1.0-0.00266*cos(2.0*pos(1))-0.00028*hgt/1E+3)/cos(z);
        trpw=0.002277*(1255.0/temp+0.05)*Edtrp/cos(z);
        zhd=trph+trpw;%对流层修正
    end
    
    if(pos(3)<-1000.0|pos(3)>20000.0)
        tropmapF=0;
    else
        el_t=azel(i*2);
        lat=pos(1)*R2D;
        hgt=pos(3);
        if el_t<=0
            tropmapF=0;
        else
            if lat<0
                lat_t=0.5;
            else
                lat_t=0;
            end
            y_t=(time2doy(time)-28.0)/365.25+lat_t;
            cosy=cos(2*pi*y_t);
            lat=abs(lat);
            for w=1:3
                interPC1=interpc(coef(w,:),lat);
                interPC2=interpc(coef(w+3,:),lat);
                interPC3=interpc(coef(w+6,:),lat);
                ah(w)=interPC1-interPC2*cosy;
                aw(w)=interPC3;
            end
            dm=(1/sin(el_t)-mapf(el_t,aht(1),aht(2),aht(3)))*hgt/1000;
            tropmapF=mapf(el_t,ah(1),ah(2),ah(3))+dm;
        end
    end
    r=r+tropmapF*zhd;                                                       %  修正对流层延时
        
    %设置天线相位中心修正为0
    dant=[0,0,0];
    
    %undifferenced phase/code residual for satellite
    lam(1)=CLIGHT/FREQ1; %波长
    y(2*i-1)=ObsODat_k.Obs_FreL1(i)*lam(1)-r-dant(1);%%载波相位观测量*光速*  载波相位残差  这里没用电离层模型
    y(2*i)=ObsODat_k.Obs_RangeC1(i)-r-dant(1);      %%非差分残差表达式       伪距残差
  
end

function interPC=interpc(coef_t,lat)
i_t=fix(lat/15);
if i_t<1
    interPC=coef_t(1);
elseif i_t>4
    interPC=coef_t(5);
else
    interPC=coef_t(i_t)*(1-lat/15+i_t)+coef_t(i_t+1)*(lat/15.0-i_t);
end

function mapF=mapf(el_t,ah,bh,ch)
sinel=sin(el_t);
mapF=(1+ah/(1+bh/(1+ch)))/(sinel+(ah/(sinel+bh/(sinel+ch))));
