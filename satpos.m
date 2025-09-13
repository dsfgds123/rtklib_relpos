function [x,y,z,dts,var]= satpos(EphDatn,time)
% input:
% EphDatn             ����
% time                ����ʱ��
% *********************************
% output��
% x,y,z               ����λ��
% dts                 �����Ӳ�
% var                 ����λ�÷���
% *********************************

EphDatn.toe;%�����ο�ʱ��
tk=timediff(time,EphDatn.toe);%��    �滮ʱ��   timeһ����toeǰ�����Сʱ֮��
global CLIGHT;

mu=3.9860050E+14;   %�������������׳�GM����u
we=7.2921151467E-5;      %������ת�ٶ�
RTOL_KEPLER=1E-14;     %����Ҫ��
MAX_ITER_KEPLER=30;%�����������ֵ

a=EphDatn.sqa^2;           %���򳤰���
M=EphDatn.M0+(sqrt(mu/a^3)+EphDatn.dn)*tk;%ƽ�����       %sqrt(mu/a^3)+EphDatn.dnΪ����ƽ�����ٶ�
E=M;                                         %ƫ�����  ����
n=0;
Ek=0;
%�������ǵ�ƫ�����Ek
while(abs(E-Ek)>RTOL_KEPLER&&n<MAX_ITER_KEPLER)%��������30�������㾫��Ҫ�������
    Ek=E;
    E=E-(E-EphDatn.e*sin(E)-M)/(1-EphDatn.e*cos(E));%����
    n=n+1;
end
%���������ܳ���n>30���������������һ��while����У�n���Ŀ���ȡֵ��ֻ��30
if(n>=MAX_ITER_KEPLER)%����������ʮ��
    fprintf('kepler iteration overflow');%����Ʈ�� 
end
iuk=atan2(sqrt(1-EphDatn.e*EphDatn.e)*sin(E),cos(E)-EphDatn.e)+EphDatn.w;  %������Ǿ�
irk=a*(1-EphDatn.e*cos(E));  %����ʸ�����Ȼ�����
iik=EphDatn.i0+EphDatn.iDot*tk; %�����ǻ�����
%������+���������ڸ��������
u=EphDatn.Cuc*cos(2*iuk)+EphDatn.Cus*sin(2*iuk)+iuk;            %����������Ǿ�
r=EphDatn.Crc*cos(2*iuk)+EphDatn.Crs*sin(2*iuk)+irk;            %������ĵ����ྶ
i=EphDatn.Cic*cos(2*iuk)+EphDatn.Cis*sin(2*iuk)+iik;            %���������� 

xx=r*cos(u);                        %�ڹ�����е�λ��
yy=r*sin(u);
omg=EphDatn.omg0+(EphDatn.omg-we)*tk-we*EphDatn.toes;   %������������㾭��
x=xx*cos(omg)-yy*cos(i)*sin(omg); %���ĵع�ֱ������ϵ�е�λ�ã�x,y,z��
y=xx*sin(omg)+yy*cos(i)*cos(omg);
z=yy*sin(i);
tk=time(1)-EphDatn.tocG(1)+time(2)-EphDatn.tocG(2);  %��һ�����������Ӳ�
dts1=EphDatn.a0+EphDatn.a1*tk+EphDatn.a2*tk*tk;
dts=dts1-(2*sqrt(mu*a)*EphDatn.e*sin(E)/CLIGHT/CLIGHT);   %�����У����  ��P79

ura_value=[2.4 3.4 4.85 6.85 9.65 13.65 24 48 96 192 384 768 1536 3072 6144];
if(EphDatn.svacc<0||EphDatn.svacc>15)
    var=6144*6144;
else
    var=ura_value(EphDatn.svacc)*ura_value(EphDatn.svacc);%����λ�÷������
end

satpos=[x,y,z,dts,var];
% ��ӡ���
% x=num2str(x,'%12.6f');
% y=num2str(y,'%12.6f');
% z=num2str(z,'%12.6f');
% fprintf('���Ǻ�PRN:%2.0f',PRN); fprintf('\n');   
% fprintf('ʱ��time:%4.0f%4.0f%4.0f%4.0f%4.0f%4.1f',time); fprintf('\n');
% fprintf('���������ڵ�������ϵ�еĿռ�����Ϊ��\n');
% fprintf('X = ');fprintf(x,'\n');fprintf('m\n');
% fprintf('Y = ');fprintf(y,'\n');fprintf('m\n');
% fprintf('Z = ');fprintf(z,'\n');fprintf('m\n');
% 
return