function t=adjweek(t1,t0)

tt=timediff(t1,t0);
if(tt<-302400)
    t=timeadd(t1,604800);
elseif(tt>302400)
    t=timeadd(t1,-604800);
else
    t=t1;
end