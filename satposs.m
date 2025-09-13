function [RS,DTS,VAR]=satposs(time,tJD,ObsODatn,Snum,EphDat,R,N)%������λ�õĺ���
% input:
% time                         ����Ԫ���źŽ���ʱ��
% tJD                          �źŽ���ʱ���������
% ObsODatn                     ��n����Ԫ�Ĺ۲�����
% Snum                         ����Ԫ�۲⵽�����Ǹ���
% EphDat                       ��������
% R                            L1��α��۲�ֵ
% N                            �������ݵ����Ǹ���
% *************************************************************************
% output:
% RS                           ����λ��
% DTS                          �����Ӳ�
% VAR                          ����λ�õķ���
% *************************************************************************

global CLIGHT;%����

Code=ObsODatn.SatCode;       %����Ԫ�۲⵽�����Ǻ�

for j=1:Snum
    
    for i=1:N               %ÿ����Ԫ���۲⵽�����ǻ������ǵ��������ļ�������
        timei= timeadd(time,-R(j)/CLIGHT);%transmission time by satellite clock ���õ�������źŷ���ʱ��
        if(EphDat(i).PRN==Code(j)&&abs(tJD-EphDat(i).toc)<0.0417*2)%����ʱ������������ļ�
            %satellite clock bias by broadcast ephemeris  �ɹ㲥���������������ʱ��ƫ��
            a0= EphDat(i).a0;
            a1= EphDat(i).a1;
            a2= EphDat(i).a2;
            %toe=EphDat(i).toe;
            ta=((timei(1)-EphDat(i).tocG(1))+(timei(2)-EphDat(i).tocG(2)));%toc�Ӳ�����Ĳο�ʱ�䣬t-toc  toe �����ο�ʱ��
            for w=1:2
                ta=ta-(a0+a1*ta+a2*ta*ta);  %������������Ӳ� (a0+a1*ta+a2*ta*ta)Ϊ�����Ӳ�
            end
%             tb=ta-(a0+a1*ta+a2*ta*ta);
%             tc=tb-(a0+a1*tb+a2*tb*tb);
            dt=a0+ta*a1+ta*ta*a2;%�����Ӳ�
            
            timei=timeadd(timei,-dt);
            
            %satellite position and clock at transmission time ����λ�ú�ʱ���ڴ���ʱ��
            [rs(1),rs(2),rs(3),dts,var]= satpos(EphDat(i),timei);
            if j==1
                RS=[rs(1),rs(2),rs(3)];
            else
                RS=[RS;rs(1),rs(2),rs(3)];
            end
            DTS(j)=dts;
            VAR(j)=var;
            %ts=ts2+tm;%����Ϊ����Ư
            %[Xs2,Ys2,Zs2,dts2]= CalPos1(EphDat(i),ts2);
            %Xs=(Xs2-Xs)/tm; %Ys=(Ys2-Ys)/tm; %Zs=(Zs2-Zs)/tm;
            %dts=(dts2-dts)/tm;
            %break;
        end
    end
end
