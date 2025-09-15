function [rtk,xa,nb]=resamb_LAMBDA(rtk)
% *************************************************************************
% input:
% rtk                          ��������ṹ��
% *************************************************************************
% output:
% rtk                          ��������ṹ��
% xa                           ���ջ�״̬����
% nb                           ģ���ȸ���
% *************************************************************************

global MAXSAT;

nx=rtk.nx; %״̬��������ʼֵΪ35��                    %float״̬������
na=rtk.na; %����״̬��������ʼֵΪ3��                    %fixed״̬������
rtk.sol.ratio=0.0;             %����ģ������Ч�ԵĲ���

D_=zeros(nx);

[nb,D,SsatFix]=ddmat(rtk,D_);%SsatFix:ģ���ȱ�־��0��nodata��1��float��2��fix��3��hold
                             %�󵥲�ת˫�����
for i=1:MAXSAT
    rtk.ssat(i).fix=SsatFix(i);
end

if nb<=0      %nbģ���ȸ���
    error('no valid double-different');
end

ny=na+nb;
y=D'*rtk.x;                   %�õ�״̬�����˫��
Qy=D'*rtk.P*D;                %�õ�״̬��˫���Э�������
for i=1:nb
    for j=1:nb
        Qb(i,j)=Qy(na+i,na+j);%�õ�ģ���ȵķ���Э������QN
    end
end

for i=1:na
    for j=1:nb
        Qab(i,j)=Qy(i,na+j);%QNR
    end
end

for i=1:nb
    yna(i)=y(na+i);%˫��ģ����   ����ǰ������bias�Ĳ�ֵ��˫��ģ����
end

yna=yna';

[info,b,s]=lambda(nb,2,yna,Qb);

if info==0
    if s(1)>0

       rtk.sol.ratio=s(2)/s(1);

    else
        rtk.sol.ratio=0.0;
    end
    if rtk.sol.ratio>999.9
        rtk.sol.ratio=999.9;
    end

    %validation by popular ratio-test,ratio-test���
     if (s(1)<=0.0|s(2)/s(1)>=3.0)   %ratio��ֵ����Ϊ3       
for i=1:na
            rtk.xa(i)=rtk.x(i);   
            for j=1:na
                rtk.Pa(i,j)=rtk.P(i,j);
            end
        end
        for i=1:nb
            bias(i)=b(1,i);  %bΪ����ģ���Ƚ��������F
            y(na+i)=y(na+i)-b(1,i);
        end
        Qb=inv(Qb);%����
        for i=na+1:ny
            y1(i-na)=y(i);
        end
        y1=y1';
        db=Qb*y1;
        rtk.xa=rtk.xa-(Qab*db);  %rtk.xa  ģ���Ƚ������ջ�״̬
        QQ=Qab*Qb;
        rtk.Pa=rtk.Pa-(QQ*Qab'); %rtk.pa  ģ���Ƚ���󷽲�Э������
        
        %�ָ�����ģ����
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
%�õ���-˫��ת������D
global MAXSAT;
global D2R;

nx=rtk.nx;
na=rtk.na;%����״̬���� ��ʼֵΪ3
D=D_;

for i=1:MAXSAT
    SsatFix(i)=0;
end

for i=1:na
    D(i,i)=1.0;   %�������D
end
k=na;
nb=1;

for i=k+1:k+MAXSAT
    if rtk.x(i)==0|~(rtk.ssat(i-k).vsat)%������Ч��־���������Ǻ�  ������Ϊ0�������ߣ�
        continue;
    end
    if rtk.ssat(i-k).lock>0&~(bitand(rtk.ssat(i-k).slip,2))&rtk.ssat(i-k).azel(2)>=15*D2R %������������1����������2�����������ǣ�3�����Ǵ���15��  %15������Ӧ�Ļ���
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
  %  if rtk.ssat(j-k).lock>0&~(bitand(rtk.ssat(j-k).slip,2))&rtk.ssat(j-k).vsat&(rtk.ssat(j-k).azel(2)>=15*D2R) %�ο����ǵ�ѡ����������������
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
function [info,F,s]=lambda(n,m,a,Q)  %�����a��Ϊyna˫��ģ����
% lambda/mlambda integer least-square estimation
if(n<=0||m<=0)
    info=-1;
    error('lambda ��ֵ����');
end
L=zeros(n);
z=zeros(n,1);
E=zeros(n,m);
% LD factorization
[info,L_,D_]=LD(n,Q);%�ֽ����QΪL�����Q����

if info==0
    % lambda reduction
    [Z,L,D]=reduction(n,L_,D_);%����Z�任
    z=Z'*a;

    % mlambda search
    [info,E,s]=search(n,m,L,D,z);   
    if info==0;
        [info,F]=solve(Z,E);%EΪ���ó�����������ģ�����������Ž�ʹ��Ž⣬ZΪ�仯����
    end
    F=F';
end
% *************************************************************************
function [info,L,D]= LD(n,Q)%�� The method for integer ambiguity implementation aspects �еĵ������н���
%�ú�������״̬����Э����������ʽ�ֽ�  ��Q������зֽ�
info=0;
A=zeros(n);
D=zeros(n,1);
A=Q;
for i=n:-1:1
    D(i)=A(i,i);
    if D(i)<=0.0   %Ϊʲô�ж�������С��0
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
function [Z,L,D]=reduction (n,L_,D_)%�����
%����Z�任
L=L_;
D=D_;
j=n-2+1;%7
k=n-2+1;%7
Z=eye(n);%��λ��
while j>=1
    if j<=k
        for i=j+1:n
            [Z,L]=gauss(n,L,Z,i,j);
        end
    end
    
    del=D(j)+L(j+1,j)*L(j+1,j)*D(j+1);%delΪD��j+1���������ֵ
    
    if (del+1e-6)<D(j+1)%���������ֵ��ԭ����ֵ���Ƚ�
        [Z,L,D]=perm(n,L,D,j,del,Z);
        k=j;
        j=n-1;
    else
        j=j-1;
    end
end
% *************************************************************************
function [Zi,L]=gauss(n,L,Zi,i,j)  %GPS���������ݴ���P175
%���и�˹�任
mu=round(L(i,j));%�ж��Ƿ�С�ڵ���0.5
if(mu~=0)
    for k=i:n
        L(k,j)= L(k,j)-mu*L(k,i);
    end
    for k=1:n
        Zi(k,j)=Zi(k,j)-mu*Zi(k,i);%Z��λ�任����ZiΪ���������ı仯����
    end
end
% *************************************************************************
function [Z,L,D]=perm(n,L,D,j,del,Z)
%�����û�   ��ZΪ�н���
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
for k=1:n    %��Z�����໥�û�
    h=Z(k,j);
    Z(k,j)=Z(k,j+1);
    Z(k,j+1)=h;
end
% *************************************************************************
function [info,zn,s]=search(n,m,L,D,zs)
%��������ģ���ȵ�����
imax=1;
nn=1;%��ѡ������z�ĸ���
maxdist=1e99;%������
S=zeros(n,n);
dist=zeros(n,1);%�ۺ�
zb=zeros(n,1);
z=zeros(n,1);
step=zeros(n,1);

zn=zeros(n,m);%��ѡ������z
s=zeros(2,1);

k=n;
dist(k)=0.0;
zb(k)=zs(k);%����
z(k)=round(zb(k));%��
y=zb(k)-z(k);%��ֵ�����ں����ʾ����
if y<=0.0
    step(k)=-1.0;
else
    step(k)= 1.0;  
end

for c=1:10000
    %����ÿ����Ҫ��������k=1��֮��Ż�����else��306����ÿ��������k~=1�������newdist���и��£���k=1ʱ���Ѿ�����������ۼƺ͵ĵ��ӡ�
    newdist=dist(k)+y*y*1/D(k);%����
    if newdist<maxdist
        if k~=1
            k=k-1;
            dist(k)=newdist;
            for i=1:k
                S(k,i)=S(k+1,i)+(z(k+1)-zb(k+1))*L(k+1,i);%Ϊ�˵�����Zb
            end
            zb(k)=zs(k)+S(k,k);
            z(k)=round(zb(k));
            y=zb(k)-z(k);
            if(y<=0)   %��ʼ����һ�׶�step��ֵ�������ֵ�����ʲ��ܻع���һ�׶Ρ�
                step(k)=-1.0;
            else
                step(k)=1.0;
            end
        else%�����д˺�����ʱ���Ѿ�����������ۼƺͣ���ʼ��z1ֵ���и��ģ���������ѡ�⣬
            if nn<=m %�������Ž�ʹ��Ž��Ӧ��2���Ⱥ���ߵ�ֵ
                if ((nn==1)|(newdist>s(imax)))
                    imax=nn;
                end
                for i=1:n
                    zn(i,nn)=z(i);
                end
                s(nn)=newdist;%ratio����еķ��Ӻͷ�ĸ���������ƽ�����Է���
                nn=nn+1;
            else%��������ѡ����ۼƺ�û�г���maxdist���ͽ��������ۼƺͽ��бȽϣ������ʱ���ۼƺ�С��֮ǰ���ۼƺͣ���ֵ������newdist��
                %���򣬸������maxdistΪ�ڶ����ѡ����ۼƺ͡�
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
            z(1)=z(1)+step(1);%��һ���Ҿ��ÿ�����΢���´���
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
for i=1:m-1  %������ѡ������Zֵ���õ���Ӧ�Ŀ���ֵ����Ƚ�ѡ�����Ž⣬��s(1)Ϊ���Ž�
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
  info=-1; %���λ����˼��ѭ����һ���Ҳû�еó����
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
%�ָ�����ģ����
global MAXSAT;
nf=1;
nv=1;
for i=1:rtk.nx
    xa(i)=rtk.x(i);
end
for i=1:rtk.na
    xa(i)=rtk.xa(i);%��ģ���Ƚ����Ľ��ջ�״̬������ԭ����״̬
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