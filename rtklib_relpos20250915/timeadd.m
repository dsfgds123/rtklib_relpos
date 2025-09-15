function t=timeadd(time, sec)  %time当前时间化为秒和不足秒

tm(2)=time(2)+sec;
tt=floor(tm(2));%向下取整
t(1)=time(1)+tt;
t(2)=tm(2)-tt;%进位