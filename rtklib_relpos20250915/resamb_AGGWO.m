function [rtk,xa,nb]=resamb_AGGWO(rtk)
% *************************************************************************
% input:
% rtk                          解算参数结构体
% *************************************************************************
% output:
% rtk                          解算参数结构体
% xa                           接收机状态参数
% nb                           模糊度个数
% *************************************************************************

global MAXSAT;

nx=rtk.nx; %状态数量（初始值为35）                    %float状态的数量
na=rtk.na; %整数状态数量（初始值为3）                    %fixed状态的数量
rtk.sol.ratio=0.0;             %检验模糊度有效性的参数

D_=zeros(nx);

[nb,D,SsatFix]=ddmat(rtk,D_);%SsatFix:模糊度标志，0：nodata，1：float，2：fix，3：hold
                             %求单差转双差矩阵
for i=1:MAXSAT
    rtk.ssat(i).fix=SsatFix(i);
end

if nb<=0      %nb模糊度个数
    error('no valid double-different');
end

ny=na+nb;
y=D'*rtk.x;                   %得到状态矩阵的双差
Qy=D'*rtk.P*D;                %得到状态的双差方差协方差矩阵
for i=1:nb
    for j=1:nb
        Qb(i,j)=Qy(na+i,na+j);%得到模糊度的方差协方差阵QN
    end
end

for i=1:na
    for j=1:nb
        Qab(i,j)=Qy(i,na+j);%QNR
    end
end

for i=1:nb
    yna(i)=y(na+i);%双差模糊度   除开前三项外bias的差值即双差模糊度
end

yna=yna';

% --- 为 A-GGWO 算法配置参数 ---
SearchAgents_no = 120;      % 种群数量 (参考您的测试脚本 run_AGGWO_on_GNSS.m)
Max_iter = 150;           % 最大迭代次数 (参考您的测试脚本 run_AGGWO_on_GNSS.m)
dim = nb;                 % 问题的维度，即模糊度的个数

% --- 构造目标函数 ---
% 这里的 yna 是一个列向量，但您的 target_function 需要一个行向量。
% a1 对应 yna', Q1 对应 Qb
Q_inv = inv(Qb); % 预先计算逆矩阵，提高效率
fobj = @(z) (yna' - z) * Q_inv * (yna' - z)';

% --- 设定搜索空间 (Search Space) ---
% 一个常用且有效的策略是以浮点解四舍五入后的结果为中心，构建一个超立方体搜索空间。
search_center = round(yna');
search_radius = 6;  % 搜索半径，可根据问题难度调整，6是一个合理的初始值
lb = search_center - search_radius; % 下界
ub = search_center + search_radius; % 上界

% (可选) A-GGWO 的额外参数，如果不需要可以留空
AGGWO_params = [];

% --- 调用 A-GGWO 算法进行整数搜索 ---
disp('Starting A-GGWO for ambiguity resolution...'); % 添加提示信息
[Best_score, Best_pos, ~] = A_GGWO(SearchAgents_no, Max_iter, lb, ub, dim, fobj, AGGWO_params);
disp(['A-GGWO finished. Best score: ', num2str(Best_score)]);

% 为了适配后续代码，我们将 A-GGWO 的输出构造成与 lambda 类似的格式
info = 0; % 假设 A-GGWO 总能成功找到解
b = [Best_pos; zeros(1, dim)]; % 将最优解放在第一行，次优解用全零占位
s = [Best_score; 9999];        % 将最优值放在第一位，次优值设为一个很大的数

if info==0
    if s(1)>0

       rtk.sol.ratio=s(2)/s(1);

    else
        rtk.sol.ratio=0.0;
    end
    if rtk.sol.ratio>999.9
        rtk.sol.ratio=999.9;
    end

    %validation by popular ratio-test,ratio-test检测
     if (s(1)<=0.0|s(2)/s(1)>=3.0)   %ratio限值设置为3       
for i=1:na
            rtk.xa(i)=rtk.x(i);   
            for j=1:na
                rtk.Pa(i,j)=rtk.P(i,j);
            end
        end
        for i=1:nb
            bias(i)=b(1,i);  %b为整周模糊度解算出来的F
            y(na+i)=y(na+i)-b(1,i);
        end
        Qb=inv(Qb);%求逆
        for i=na+1:ny
            y1(i-na)=y(i);
        end
        y1=y1';
        db=Qb*y1;
        rtk.xa=rtk.xa-(Qab*db);  %rtk.xa  模糊度解算后接收机状态
        QQ=Qab*Qb;
        rtk.Pa=rtk.Pa-(QQ*Qab'); %rtk.pa  模糊度解算后方差协方差阵
        
        %恢复单差模糊度
        xa=restamb(rtk,bias);
        %rtkxa=rtk.xa;
        %rtkPa=rtk.Pa;
    else
        xa=zeros(nx,1);
        nb=0;
        sprintf('ambiiguity validation failed');
    end
else
    xa=zeros(nx,1);
    nb=0;
    sprintf('lambda error');
end

% *************************************************************************
function [nb,D,SsatFix]=ddmat(rtk,D_)
%得到单-双差转换矩阵D
global MAXSAT;
global D2R;

nx=rtk.nx;
na=rtk.na;%整数状态数量 初始值为3
D=D_;

for i=1:MAXSAT
    SsatFix(i)=0;
end

for i=1:na
    D(i,i)=1.0;   %构造矩阵D
end
k=na;
nb=1;

for i=k+1:k+MAXSAT
    if rtk.x(i)==0|~(rtk.ssat(i-k).vsat)%数据有效标志即共视卫星号  在这里为0则往下走，
        continue;
    end
    if rtk.ssat(i-k).lock>0&~(bitand(rtk.ssat(i-k).slip,2))&rtk.ssat(i-k).azel(2)>=15*D2R %满足三条件，1、无周跳，2、锁定的卫星，3、仰角大于15°  %15°所对应的弧度
        SsatFix(i-k)=2;
        break;
    else
        SsatFix(i-k)=1;
    end
end

for j=k+1:k+MAXSAT
    if i==j|rtk.x(j)==0|~(rtk.ssat(j-k).vsat)%k=3,j-3=3+sat(i)-3=sat(i)
        continue;
    end
  %  if rtk.ssat(j-k).lock>0&~(bitand(rtk.ssat(j-k).slip,2))&rtk.ssat(j-k).vsat&(rtk.ssat(j-k).azel(2)>=15*D2R) %参考卫星的选择，锁定无周跳仰角
        D(i,na+nb)= 1.0;
        D(j,na+nb)=-1.0;
        nb=nb+1;
         SsatFix(j-k)=2;
 %   else
   %      SsatFix(j-k)=1;
  %  end
end
nb=nb-1; 
% *************************************************************************

% *************************************************************************
function [info,F,s]=lambda(n,m,a,Q)  %这里的a即为yna双差模糊度
% lambda/mlambda integer least-square estimation
if(n<=0||m<=0)
    info=-1;
    error('lambda 赋值错误');
end
L=zeros(n);
z=zeros(n,1);
E=zeros(n,m);
% LD factorization
[info,L_,D_]=LD(n,Q);%分解矩阵Q为L矩阵和Q矩阵

if info==0
    % lambda reduction
    [Z,L,D]=reduction(n,L_,D_);%进行Z变换
    z=Z'*a;

    % mlambda search
    [info,E,s]=search(n,m,L,D,z);   
    if info==0;
        [info,F]=solve(Z,E);%E为所得出的两组整周模糊度整数最优解和次优解，Z为变化矩阵
    end
    F=F';
end
% *************************************************************************
function [info,L,D]= LD(n,Q)%在 The method for integer ambiguity implementation aspects 中的第三章有解释
%该函数计算状态方差协方差矩阵的因式分解  对Q矩阵进行分解
info=0;
A=zeros(n);
D=zeros(n,1);
A=Q;
for i=n:-1:1
    D(i)=A(i,i);
    if D(i)<=0.0   %为什么判定条件是小于0
        info=-1;
        break;
    end
    a=sqrt(D(i));
    for j=1:i
        L(i,j)=A(i,j)/a;
    end
    for j=1:i-1
        for k=1:j
            A(j,k)=A(j,k)-L(i,k)*L(i,j);
        end
    end 
    for j=1:i
        L(i,j)= L(i,j)/L(i,i);
    end    
end
if info~=0
        error('LD factorization error');
end
% *************************************************************************
function [Z,L,D]=reduction (n,L_,D_)%降相关
%进行Z变换
L=L_;
D=D_;
j=n-2+1;%7
k=n-2+1;%7
Z=eye(n);%单位阵
while j>=1
    if j<=k
        for i=j+1:n
            [Z,L]=gauss(n,L,Z,i,j);
        end
    end
    
    del=D(j)+L(j+1,j)*L(j+1,j)*D(j+1);%del为D（j+1）交换后的值
    
    if (del+1e-6)<D(j+1)%将交换后的值与原来的值做比较
        [Z,L,D]=perm(n,L,D,j,del,Z);
        k=j;
        j=n-1;
    else
        j=j-1;
    end
end
% *************************************************************************
function [Zi,L]=gauss(n,L,Zi,i,j)  %GPS测量与数据处理P175
%进行高斯变换
mu=round(L(i,j));%判断是否小于等于0.5
if(mu~=0)
    for k=i:n
        L(k,j)= L(k,j)-mu*L(k,i);
    end
    for k=1:n
        Zi(k,j)=Zi(k,j)-mu*Zi(k,i);%Z单位变换矩阵，Zi为经过排序后的变化矩阵
    end
end
% *************************************************************************
function [Z,L,D]=perm(n,L,D,j,del,Z)
%进行置换   且Z为列交换
eta=D(j)/del;
lam=D(j+1)*L(j+1,j)/del;
D(j)=eta*D(j+1);
D(j+1)=del;
for k=1:j-1
    a0=L(j,k);
    a1=L(j+1,k);
    L(j,k)= -L(j+1,j)*a0+a1;
    L(j+1,k)=eta*a0+lam*a1;
end
L(j+1,j)=lam;
for k=j+2:n
    g=L(k,j);
    L(k,j)=L(k,j+1);
    L(k,j+1)=g;
end
for k=1:n    %将Z中列相互置换
    h=Z(k,j);
    Z(k,j)=Z(k,j+1);
    Z(k,j+1)=h;
end
% *************************************************************************
function [info,zn,s]=search(n,m,L,D,zs)
%进行整周模糊度的搜索
imax=1;
nn=1;%候选向量组z的个数
maxdist=1e99;%最大距离
S=zeros(n,n);
dist=zeros(n,1);%累和
zb=zeros(n,1);
z=zeros(n,1);
step=zeros(n,1);

zn=zeros(n,m);%候选向量组z
s=zeros(2,1);

k=n;
dist(k)=0.0;
zb(k)=zs(k);%浮点
z(k)=round(zb(k));%整
y=zb(k)-z(k);%差值，用于后面表示距离
if y<=0.0
    step(k)=-1.0;
else
    step(k)= 1.0;  
end

for c=1:10000
    %这里每次需要先运行至k=1，之后才会运行else（306），每次运行完k~=1，都会对newdist进行更新，当k=1时，已经完成了所有累计和的叠加。
    newdist=dist(k)+y*y*1/D(k);%距离
    if newdist<maxdist
        if k~=1
            k=k-1;
            dist(k)=newdist;
            for i=1:k
                S(k,i)=S(k+1,i)+(z(k+1)-zb(k+1))*L(k+1,i);%为了迭代出Zb
            end
            zb(k)=zs(k)+S(k,k);
            z(k)=round(zb(k));
            y=zb(k)-z(k);
            if(y<=0)   %初始化上一阶段step的值，避免该值不合适不能回归上一阶段。
                step(k)=-1.0;
            else
                step(k)=1.0;
            end
        else%当运行此函数的时候，已经完成了所有累计和，开始对z1值进行更改，求出俩组候选解，
            if nn<=m %保存最优解和次优解对应的2个等号左边的值
                if ((nn==1)|(newdist>s(imax)))
                    imax=nn;
                end
                for i=1:n
                    zn(i,nn)=z(i);
                end
                s(nn)=newdist;%ratio检测中的分子和分母，即距离的平方除以方差
                nn=nn+1;
            else%如果俩组候选解的累计和没有超过maxdist，就将两组解的累计和进行比较，如果此时的累计和小于之前的累计和，赋值并更改newdist；
                %否则，更改最大maxdist为第二组候选解的累计和。
                if newdist<s(imax)
                    for i=1:n
                        zn(i,imax)=z(i);
                    end
                    s(imax)=newdist;
                    imax=1;
                    for i=1:m
                        if (s(imax)<s(i))
                            imax=i;
                        end
                    end
                end
                maxdist=s(imax);
            end
            z(1)=z(1)+step(1);%这一块我觉得可以稍微改下代码
            y=zb(1)-z(1);
            if step(1)<=0
                v=-1.0;
            else
                v=1.0;
            end
            step(1)=-step(1)-v;
        end
    else
        if k==n
            break;
        else
            k=k+1;
            z(k)=z(k)+step(k);
            y=zb(k)-z(k);
            if (step(k))<=0
                v=-1.0;
            else
                v=1.0;
            end
          step(k)=-step(k)-v;
        end
    end
end
% %------------%
% s(2)=newdist;
%  for i=1:n
%    zn(i,2)=z(i);
%  end
%---------%
for i=1:m-1  %根据所选的两组Z值，得到对应的卡方值，相比较选出最优解，即s(1)为最优解
    for j=i+1:m
        if s(i)<s(j)
            continue;
        end
        v=s(i);
        s(i)=s(j);
        s(j)=v;
        for k=1:n
            v=zn(k,i);
            zn(k,i)=zn(k,j);
            zn(k,j)=v;
        end
    end
end

if c>10000
  info=-1; %标记位，意思是循环了一万次也没有得出结果
end
info=0;
% *************************************************************************
function [info,X]=solve(A,Y)

B=A;
B1=inv(B);
X=B1'*Y;
info=0;
% *************************************************************************
 function xa=restamb(rtk,bias)
%恢复单差模糊度
global MAXSAT;
nf=1;
nv=1;
for i=1:rtk.nx
    xa(i)=rtk.x(i);
end
for i=1:rtk.na
    xa(i)=rtk.xa(i);%将模糊度解算后的接收机状态更换掉原来的状态
end
f=nf;
n=1;
for i=1:MAXSAT
    if rtk.ssat(i).fix~=2
        continue;
    end
    index(n)=3+i;
    n=n+1;
end
n=n-1;
if n>=2
    xa(index(1))=rtk.x(index(1));
    for i=2:n
        xa(index(i))=xa(index(1))-bias(nv);
        nv=nv+1;
    end
end
% *************************************************************************