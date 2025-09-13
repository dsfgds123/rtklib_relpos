function [nv,rtk,v,H,R,vflg]=ddres(rtk,dt,x,sat,y_r,e_r,azel_r,y_b,ir,ib,ns,Hmode)
% *************************************************************************
% input:
% rtk                          ��������ṹ��
% dt                           ��վ�ƶ�վ������������ʱ���
% x                            �ƶ�վ���ջ�λ��
% sat                          ��ʽ���ǵ����Ǻ�
% y_r                          �ƶ�վ���ջ��ǲ�ֲв����α����ز���λ��
% e_r                          �ƶ�վ���ջ�������֮�������
% azel_r                       �ƶ�վ���ջ�������֮��ĸ߶ȽǺ�����
% y_b                          ��վ���ջ��ǲ�ֲв�
% ir                           �������ǵ��ƶ�վ���ջ����Ǻ�����
% ib                           �������ǵĻ�վ���ջ����Ǻ�����
% ns                           �������Ǹ���
% Hmode                        ��������ϵ������H
% *************************************************************************
% output:
% nv                           ˫��в����
% rtk                          ��������ṹ��
% v                            ˫��в�����
% H                            ϵ������H
% R                            ˫��������ķ���Э����
% vflg                         ��־λ
% *************************************************************************

global CLIGHT;
global MAXSAT;
global FREQ1;

H=zeros(ns*2,rtk.nx);%%����������Ŀ*2��״̬����������һ����18��35���������

didxi=0;
didxj=0;
nv=1;
b=1;     %%����
nb=[0,0];%%����

[bl,dr]=baseline(x,rtk.rb);%x�ƶ�վλ�ú�rb��վλ��  bl���߾���   dr����ʸ��

posu=ecef2pos(x);
posr=ecef2pos(rtk.rb);%תΪ��γ��

for i=1:MAXSAT
    rtk.ssat(i).resp=0;          %rtk.ssatΪ����״̬�ṹ�壬rtk.ssat.respα��˫��вrtk.ssat.resc �ز���λ˫��в�
    rtk.ssat(i).resc=0;
end

for f=1:2
    %search reference satellite with highest elevation
    %��������ߵ����ǵ����ο����Ǽ�������ߵĹ�������
    i=-1;
    for j=1:ns
        if i<0|azel_r(2*ir(j))>=azel_r(2*ir(i))
            i=j;
        end
    end    %ȷ��i��ֵ��ȷ��������ߵĹ������ǣ��������õ��������������i=3��
    if i<0
        continue;
    end
    
    %make double difference
    for j=1:ns
        if i==j
            continue;
        end

        lami=CLIGHT/FREQ1;  %%����
        lamj=CLIGHT/FREQ1;
        
        %double-differenced residual%   f=1 �õ��ز���λ˫��вf=2�õ�α��˫��в�
        v(nv)=(y_r(f+(ir(i)-1)*2)-y_b(f+(ib(i)-1)*2))...   %����������վ��������Ǽ��
            -(y_r(f+(ir(j)-1)*2)-y_b(f+(ib(j)-1)*2)); 
        
        %partial derivatives by rover position �ƶ�վλ��ƫ��
        if Hmode==1
            for k=1:3
                H(nv,k)=-e_r(ir(i),k)+e_r(ir(j),k);%nvΪ�в�������H�����ĵ�����ϸ�ļ���
            end
        end
        
        %double-differenced phase-bias term ˫����λ�в����ģ���ȵ����
        if f==1
            v(nv)=v(nv)-lami*x(3+sat(i))+lami*x(3+sat(j));%  ����ֵ �� Ԥ��ֵ   ���������Ǹ���ֵ���Ҿ����Ǹ�˫��в�����
            if Hmode==1
                H(nv,3+sat(i))= lami;
                H(nv,3+sat(j))=-lamj;
            end
        end
        
        if f==1
            rtk.ssat(sat(j)).resc=v(nv);% resc�ز���λ�в�
        else
            rtk.ssat(sat(j)).resp=v(nv);% respα��в� ����RTK�ṹ��
        end
        
        %test innovation
        % reject threshold of innovation, =30
        if abs(v(nv))>30
            if f==1
                rtk.ssat(sat(i)).rejc=rtk.ssat(sat(i)).rejc+1;%�ų������ݵĸ���
                rtk.ssat(sat(j)).rejc=rtk.ssat(sat(j)).rejc+1;
            end
%             continue;
        end
        
        %single-differenced measurement error variances �����������
        Ri(nv)=vareer(azel_r(ir(i)*2),bl,dt,f);%dt��վ���ƶ�վ�������ݵ�ʱ���  ����
        Rj(nv)=vareer(azel_r(ir(j)*2),bl,dt,f);
        
        %set valid data flags ������Ч���ݱ�־
        if f==1
            rtk.ssat(sat(i)).vsat=1;
            rtk.ssat(sat(j)).vsat=1;
        end
        
        if f==1
            w=0;
        else
            w=1;
        end
        vflg(nv)=bitor(bitor(bitor(bitshift(sat(i),16),bitshift(sat(j),8)),bitshift(w,4)),0);
        nv=nv+1;
        nb(b)=nb(b)+1;
    end
    b=b+1;
end
b=b-1;
nv=nv-1;

%double-differenced measurement error covariance
R=ddcov(nb,b,Ri,Rj,nv);


function [bl_,dr_]=baseline(rr,rb)
for i=1:3
    dr_(i)=rr(i)-rb(i);%����
end
bl_=norm(dr_);%�����������ֵ��վ���ƶ�վ��ֱ�߾��루���ߣ�


function Ri_nv=vareer(el,bl,dt,f)  %�������Э�����
CLIGHT=299792458.0;
c=0*bl/1e4;%c=0
d=CLIGHT*5e-12*dt;
fact=1;

if f==2
    fact=100;
end
if fact<=0
    fact=100;
end

%opt->err[1]=opt->err[2]=0.003
a=fact*0.003;
b=fact*0.003;

Ri_nv=2*(a*a+b*b/sin(el)/sin(el)+c*c)+d*d;%�������Э������

function R_=ddcov(nb_,b_,Ri_,Rj_,nv_)  %˫�����Э�����
k=0;
R_=zeros(nv_);

for b=1:b_
    for i=1:nb_(b_)
        for j=1:nb_(b_)
            if i==j
                R_(k+i,k+j)=Ri_(k+i)+Rj_(k+i);
            else
                R_(k+i,k+j)=Ri_(k+i);
            end
        end
    end
    k=k+nb_(b);
end
