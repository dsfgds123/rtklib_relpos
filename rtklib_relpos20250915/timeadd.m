function t=timeadd(time, sec)  %time��ǰʱ�仯Ϊ��Ͳ�����

tm(2)=time(2)+sec;
tt=floor(tm(2));%����ȡ��
t(1)=time(1)+tt;
t(2)=tm(2)-tt;%��λ