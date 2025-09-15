function [SsatFix,rtkX,rtkP]=holdamb(rtk,xa)

global D2R;
global MAXSAT;

nb=rtk.nx-rtk.na;%32
nv=1;
H=zeros(nb,rtk.nx);
n=1;

for i=1:MAXSAT
    SsatFix(i)=rtk.ssat(i).fix;
end

for i=1:MAXSAT
    if rtk.ssat(i).fix~=2|rtk.ssat(i).azel<15*D2R     
        continue;
    end
    index(n)=3+i;%卫星号加三
    SsatFix(i)=3;%从fix改为hold
    n=n+1;
end

n=n-1;

for i=2:n
    v(nv)=(xa(index(1))-xa(index(i)))-(rtk.x(index(1))-rtk.x(index(i)));%选第一颗卫星为参考卫星
    H(index(1),nv)= 1;
    H(index(i),nv)=-1;
    nv=nv+1;
end
nv=nv-1;

if nv>1
    R=zeros(nv);
    for i=1:nv
        R(i,i)=0.001;
    end
    
    [rtkX,rtkP]=EKF_filter(rtk.x,rtk.P,H,v,R,rtk.nx,nv,1);
    
end
