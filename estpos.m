function sol=estpos(t,solkrr,rs,dts,Snum,R,Code,vare,EphDat)
% *************************************************************************
% input:
% t                            �źŽ���ʱ��
% solkrr                       ��һ��Ԫ�Ľ��ջ�λ��
% rs                           ����λ��
% dts                          �����Ӳ�
% Snum                         ������Ŀ
% R                            α��۲�ֵ
% Code                         ���Ǻ�����
% vare                         ����λ�õķ���
% EphDat                       ��������
% *************************************************************************
% output:
% sol                          ���㶨λ�Ľ�
% *************************************************************************

global CLIGHT;

 X=[solkrr(1);solkrr(2);solkrr(3);0];

 for i=1:10
     [v,var,nv,H]=PNTrescode(X,t,rs,dts,vare,Snum,R,Code,EphDat);%�����������ǵĹ���
%      if(nv<4)
%          continue;
%      end
     %weight by variance
     for j=1:nv
         sig=sqrt(var(j));%��Ȩ��С����
         v(j)=v(j)/sig;
         H(j,:)=H(j,:)/sig;
     end
     
     %least square estimation
     v=v';
     dx=(H'*H)^-1*H'*v;
     
     X=X+dx;
     
     if(norm(dx)<1E-4)
         sol.time=timeadd(t,-X(4,1)/CLIGHT);
         sol.dtr=X(4);
         sol.rr=[X(1,1),X(2,1),X(3,1)];
         sol.age=0;
         sol.ratio=0;
         break;
     end
 end
 