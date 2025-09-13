function t=TimetoJD(Y,M,D,h,f,s)
if(Y>=80)
    Y=Y+1900;
else
    Y=Y+2000;
end
if M<=2
    Y=Y-1;
    M=M+12;
end
JD=fix(365.25*Y)+fix(30.6001*(M+1))+D+h/24+f/1440+s/24/3600+1720981.5;
t=JD-2444244.5;
%JD=Julian Date（儒略日），儒略历以天数计量，在GPS坐标转换的计算中更方便。
%GPS标准时间JD=2444244.5（1980年1月6日0时），GPS日和周的计算公式是N=mod（INT（JD+1.5），7）,W=INT（（JD-2444244.5）/7）