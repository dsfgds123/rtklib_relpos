function [xp_,Pp_]=EKF_filter(x,P,H,v,R,n,m,HOLD)
% input
% x :   ���ջ�λ��
% P :   rtk.p���������
% H :   ϵ������H��Ҳ�������ϵ�C
% v :   ˫��в�������Ҳ����˫������ģ���ȼ�����
% R :   ˫��������ķ���Э�������
% n :   rtk.nx״̬����
% m :   ��������ṹ��
% HOLD :0
% ***************************
% output
% xp_ : ���ջ�λ��
% Pp_ : �˲����P
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



%�Ѵ�����Ϣ��״̬�����䷽��Э�������ٳ��������������˲�
for i=1:k
    x_(i)=x(ix(i));
    for j=1:k
        P_(i,j)=P(ix(i),ix(j));%��i=jʱ���Ż���ֵ���룬���඼Ϊ0
    end
end

if HOLD==0  %hold=0
    for j=1:m
        for i=1:k
            H_(j,i)=H(j,ix(i));%��3+sat(i)����ֵ�������н���
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

v=v';%˫������ģ����ת��

I=eye(k);%%��λ����
x_=x_';%ת��

K=P_*H_'*(H_*P_*H_'+R)^-1;
x_=x_+K*v;%xp=x+K*v;
Pp=(I-K*H_)*P_;

%�ٽ��˲������ֵ�ش���ԭ����״̬���;���P��
for i=1:k
    xp_(ix(i))=x_(i);
    for j=1:k
        Pp_(ix(i),ix(j))=Pp(i,j);%����ʼ�ĸ�ʽ
    end
end
