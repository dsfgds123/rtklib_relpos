function enu=eceftoenu(pos,r)
%经纬高转东北天
E(1,1)=-sin(pos(2));
E(1,2)=cos(pos(2));
E(1,3)=0;

E(2,1)=-sin(pos(1))*cos(pos(2));
E(2,2)=-sin(pos(1))*sin(pos(2));
E(2,3)=cos(pos(1));

E(3,1)=cos(pos(1))*cos(pos(2));
E(3,2)=cos(pos(1))*sin(pos(2));
E(3,3)=sin(pos(1));

enu=E*r';