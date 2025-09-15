function t=time2doy(gt)
ep=time2epoch(gt);
if ep(1)>=2000
    ep(1)=ep(1)-2000;
end
ep(2)=1;
ep(3)=1;
ep(4)=0;
ep(5)=0;
ep(6)=0;
t=timediff(gt,epoch2time(ep(1),ep(2),ep(3),ep(4),ep(5),ep(6)))/86400+1;