function [xp_,Pp_]=EKF_filter(x,P,H,v,R,n,m,HOLD)
% input
% x :   接收机位置
% P :   rtk.p均方误差阵
% H :   系数矩阵H，也就是书上的C
% v :   双差残差向量，也就是双差整周模糊度加噪声
% R :   双差测量误差的方差协方差矩阵
% n :   rtk.nx状态数量
% m :   解算参数结构体
% HOLD :0
% ***************************
% output
% xp_ : 接收机位置
% Pp_ : 滤波后的P
% ***************************
xp_=x;
Pp_=P;

k=1;
for i=1:n
    if x(i)~=0&P(i,i)>-1 %3+sat(i)
        ix(k)=i;
        k=k+1;
    end
end
k=k-1;



%把带有信息的状态量和其方差协方差矩阵抠出来，作卡尔曼滤波
for i=1:k
    x_(i)=x(ix(i));
    for j=1:k
        P_(i,j)=P(ix(i),ix(j));%当i=j时，才会有值赋入，其余都为0
    end
end

if HOLD==0  %hold=0
    for j=1:m
        for i=1:k
            H_(j,i)=H(j,ix(i));%把3+sat(i)赋的值重新排列紧凑
        end
    end
elseif HOLD==1
    for i=1:k
        for j=1:m
            H_(i,j)=H(ix(i),j);
        end
    end
    H_=H_';
end

v=v';%双差整周模糊度转置

I=eye(k);%%单位矩阵
x_=x_';%转置

K=P_*H_'*(H_*P_*H_'+R)^-1;
x_=x_+K*v;%xp=x+K*v;
Pp=(I-K*H_)*P_;

%再将滤波后的数值回代到原来的状态量和矩阵P中
for i=1:k
    xp_(ix(i))=x_(i);
    for j=1:k
        Pp_(ix(i),ix(j))=Pp(i,j);%变回最开始的格式
    end
end
