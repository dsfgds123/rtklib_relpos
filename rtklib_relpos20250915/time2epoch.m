function ep=time2epoch(gt)

mday=[31,28,31,30,31,30,31,31,30,31,30,31,31,28,31,30,31,30,31,31,30,31,30,31,31,29,31,30,31,30,31,31,30,31,30,31,31,28,31,30,31,30,31,31,30,31,30,31];

days=fix(gt(1)/86400);
sec=fix(gt(1)-days*86400);
day=mod(days,1461);
for mon=1:48
    if day>=mday(mon)
        day=day-mday(mon);
    else
        break;
    end
end
ep(1)=1970+fix(days/1461)*4+fix(mon/12);
ep(2)=mod(mon,12);
ep(3)=day+1;
ep(4)=fix(sec/3600);
ep(5)=fix(mod(sec,3600)/60);
ep(6)=mod(sec,60)+gt(2);