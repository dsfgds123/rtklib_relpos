% EphDatn.svacc=[1];
% ura_value=[2.4 3.4 4.85 6.85 9.65 13.65 24 48 96 192 384 768 1536 3072 6144];
% if(EphDatn.svacc<0||EphDatn.svacc>15)
%     var=6144*6144;
% else
%     var=ura_value(EphDatn.svacc)*ura_value(EphDatn.svacc);%卫星位置方差矩阵
% end
clear;clc;
c=299792.458; %光速
%通过星历参数解算到所在地可见卫星的坐标位置
sat13=[-7134.529244 16113.648836 23709.205570];
sat22=[-22383.700040 18533.233168 5307.245613];
sat23=[-5384.901317 28971.622323 2079.796362];
sat14=[637.466571 28016.053841 9347.297933];
sat12=[-11568.199533 -3328.511543 26977.312423];
sat21=[-28908.916747 -577.061760 6051.375658];
sat5=[-1205.651181 28296.890128 -8397.025036];
sat4=[16456.527324 12347.282494 21199.173063];
sat=[sat13;sat22;sat23; sat14; sat12; sat21; sat5; sat4];%多卫星位置矩阵
%所在地实际大地坐标，用来与定位结果作比较
nanjing=[-2604.298533 4743.297217 3364.978513];
%理想伪距测量值
r0 = zeros(8,1);
for i = 1:8
    r0(i,1)=norm(sat(i,:)-nanjing);
end
%加入零均值，方差为80的随机高斯分布，模拟含有误差的伪距r
r = r0 + sqrt(80)*randn(size(r0))*1e-3;
%--------------------------------------------------------------------------
%Newton迭代法
%设定迭代初值，若无法估计则全部假设为0
xyzt=[0 0 0 0];
for i = 1:100 % 最大迭代次数设为100次
    f = dingwei_fun(xyzt,sat,r);
    df = dingwei_dfun(xyzt,sat);
    %左除，具有良好的数值稳定性，在MATLAB中此已经为最小二乘意义下的解
    %delta(xyzt)=G^(-1)*b
    delta=-df\f;
    xyzt(1)=xyzt(1)+delta(1);
    xyzt(2)=xyzt(2)+delta(2);
    xyzt(3)=xyzt(3)+delta(3);
    xyzt(4)=xyzt(4)+delta(4);
    p=norm(delta); %定位精度
    if (p<1e-10)
        break;
    end
end
disp('迭代次数为')
i
disp('用户位置为')
xyzt(1:3)
disp('用户钟差为')
xyzt(4)
%--------------------------------end---------------------------------------2、非线性方程组dingwei_fun.m
function f=dingwei_fun(xyzt,sat,r)
%多卫星定位方程函数
%xyzt(1:3)为用户位置(km)
%xyzt(4)为用户钟差(s)
%r为测量得到的伪距(km)
%sat为卫星数据
c=299792.458; %光速
[m,n]=size(sat);
for i=1:m
    f(i)=norm(sat(i,:)-xyzt(1:3))+c*xyzt(4)-r(i);
end
f=f';
end
% 3、方程组偏导数矩阵
%多颗卫星定位的偏导数矩阵
function df=dingwei_dfun(xyzt,sat)
%xyzt为用户位置及钟差
%sat为卫星数据
c=299792.458; %光速
[m,n]=size(sat);
xyz=xyzt(1:3);
for i=1:m
    for j=1:3
        df(i,j)=(xyz(j)-sat(i,j))/norm(sat(i,:)-xyz(j));
    end
end
df(:,4)=c; %线性函数ct，对t求偏导为c
end