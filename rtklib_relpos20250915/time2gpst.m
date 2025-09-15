function week=time2gpst(t)

t0=epoch2time(80,1,6,0,0,0);
sec=t(1)-t0(1);
w=fix(sec/(86400*7));

week(1)=w;
week(2)=sec-w*86400*7+t(2);
