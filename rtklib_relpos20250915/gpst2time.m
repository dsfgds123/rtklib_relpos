function t=gpst2time(week,sec)

t=epoch2time(80,1,6,0,0,0);
if(sec<-1E+9||sec>1E+9)
    sec=0.0;
end
t=[t(1)+86400*7*week+fix(sec), sec-fix(sec)];
