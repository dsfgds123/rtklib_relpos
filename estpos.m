function sol=estpos(t,solkrr,rs,dts,Snum,R,Code,vare,EphDat)
% *************************************************************************
% input:
% t                            信号接收时间
% solkrr                       上一历元的接收机位置
% rs                           卫星位置
% dts                          卫星钟差
% Snum                         卫星数目
% R                            伪距观测值
% Code                         卫星号数组
% vare                         卫星位置的方差
% EphDat                       导航数据
% *************************************************************************
% output:
% sol                          单点定位的解
% *************************************************************************

global CLIGHT;

 X=[solkrr(1);solkrr(2);solkrr(3);0];

 for i=1:10
     [v,var,nv,H]=PNTrescode(X,t,rs,dts,vare,Snum,R,Code,EphDat);%里面有求仰角的过程
%      if(nv<4)
%          continue;
%      end
     %weight by variance
     for j=1:nv
         sig=sqrt(var(j));%加权最小二乘
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
 