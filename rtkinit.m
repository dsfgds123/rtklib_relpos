function rtk=rtkinit

global MAXSAT;

rtk=struct;

rtk.nfix=0;

rtk.sol.rr=[0,0,0];
rtk.sol.qr=[0,0,0];

rtk.sol.dtr=0;%receiver clock bias
rtk.sol.stat=0;
rtk.sol.ns=0;
rtk.sol.age=0;
rtk.sol.ratio=0;
rtk.sol.time=[0,0];

rtk.rb=[0,0,0];
rtk.nx=35;
rtk.na=3;
rtk.tt=0;

rtk.x=zeros(rtk.nx,1);
rtk.P=zeros(rtk.nx);
rtk.xa=zeros(rtk.na,1);
rtk.Pa=zeros(rtk.na);

for i=1:MAXSAT
    rtk.ssat(i).resp=0;
    rtk.ssat(i).resc=0;
    rtk.ssat(i).vsat=0;
    rtk.ssat(i).fix=0;
    rtk.ssat(i).slip=0;
    rtk.ssat(i).lock=0;
    rtk.ssat(i).outc=0;
    rtk.ssat(i).rejc=0;
    rtk.ssat(i).azel=[0,0];
end
