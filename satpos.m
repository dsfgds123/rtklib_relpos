function [x,y,z,dts,var]= satpos(EphDatn,time)
% input:
% EphDatn             星历
% time                发射时刻
% *********************************
% output：
% x,y,z               卫星位置
% dts                 卫星钟差
% var                 卫星位置方差
% *********************************

EphDatn.toe;%星历参考时间
tk=timediff(time,EphDatn.toe);%减    规划时间   time一般在toe前后的两小时之间
global CLIGHT;

mu=3.9860050E+14;   %地球引力常数俗称GM或者u
we=7.2921151467E-5;      %地球自转速度
RTOL_KEPLER=1E-14;     %精度要求
MAX_ITER_KEPLER=30;%迭代次数最大值

a=EphDatn.sqa^2;           %地球长半轴
M=EphDatn.M0+(sqrt(mu/a^3)+EphDatn.dn)*tk;%平近点角       %sqrt(mu/a^3)+EphDatn.dn为卫星平均角速度
E=M;                                         %偏近点角  弧度
n=0;
Ek=0;
%计算卫星的偏近点角Ek
while(abs(E-Ek)>RTOL_KEPLER&&n<MAX_ITER_KEPLER)%迭代少于30次且满足精度要求则结束
    Ek=E;
    E=E-(E-EphDatn.e*sin(E)-M)/(1-EphDatn.e*cos(E));%迭代
    n=n+1;
end
%根本不可能出现n>30的情况，由于在上一个while语句中，n最大的可能取值就只是30
if(n>=MAX_ITER_KEPLER)%迭代超过三十次
    fprintf('kepler iteration overflow');%迭代飘了 
end
iuk=atan2(sqrt(1-EphDatn.e*EphDatn.e)*sin(E),cos(E)-EphDatn.e)+EphDatn.w;  %升交点角距
irk=a*(1-EphDatn.e*cos(E));  %卫星矢径长度基础量
iik=EphDatn.i0+EphDatn.iDot*tk; %轨道倾角基础量
%基础量+改正量等于改正后的量
u=EphDatn.Cuc*cos(2*iuk)+EphDatn.Cus*sin(2*iuk)+iuk;            %改正后的升角距
r=EphDatn.Crc*cos(2*iuk)+EphDatn.Crs*sin(2*iuk)+irk;            %改正后的地心相径
i=EphDatn.Cic*cos(2*iuk)+EphDatn.Cis*sin(2*iuk)+iik;            %改正后的倾角 

xx=r*cos(u);                        %在轨道面中的位置
yy=r*sin(u);
omg=EphDatn.omg0+(EphDatn.omg-we)*tk-we*EphDatn.toes;   %改正后的升交点经度
x=xx*cos(omg)-yy*cos(i)*sin(omg); %地心地固直角坐标系中的位置（x,y,z）
y=xx*sin(omg)+yy*cos(i)*cos(omg);
z=yy*sin(i);
tk=time(1)-EphDatn.tocG(1)+time(2)-EphDatn.tocG(2);  %又一次修正卫星钟差
dts1=EphDatn.a0+EphDatn.a1*tk+EphDatn.a2*tk*tk;
dts=dts1-(2*sqrt(mu*a)*EphDatn.e*sin(E)/CLIGHT/CLIGHT);   %相对论校正量  书P79

ura_value=[2.4 3.4 4.85 6.85 9.65 13.65 24 48 96 192 384 768 1536 3072 6144];
if(EphDatn.svacc<0||EphDatn.svacc>15)
    var=6144*6144;
else
    var=ura_value(EphDatn.svacc)*ura_value(EphDatn.svacc);%卫星位置方差矩阵
end

satpos=[x,y,z,dts,var];
% 打印结果
% x=num2str(x,'%12.6f');
% y=num2str(y,'%12.6f');
% z=num2str(z,'%12.6f');
% fprintf('卫星号PRN:%2.0f',PRN); fprintf('\n');   
% fprintf('时间time:%4.0f%4.0f%4.0f%4.0f%4.0f%4.1f',time); fprintf('\n');
% fprintf('所求卫星在地心坐标系中的空间坐标为：\n');
% fprintf('X = ');fprintf(x,'\n');fprintf('m\n');
% fprintf('Y = ');fprintf(y,'\n');fprintf('m\n');
% fprintf('Z = ');fprintf(z,'\n');fprintf('m\n');
% 
return