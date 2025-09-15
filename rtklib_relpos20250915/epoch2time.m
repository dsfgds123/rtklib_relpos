function t=epoch2time(Y,M,D,h,f,s)
if(Y>=80)
    Y=Y+1900;
else
    Y=Y+2000;
end
doy=[1 32 60 91 121 152 182 231 244 274 305 335];
days=(Y-1970)*365+fix((Y-1969)/4)+doy(M)+D-2;
if mod(Y,4)==0&(M>=3)
    days=days+1;
end
sec=s;
t=[days*86400+h*3600+f*60+sec, s-sec];
