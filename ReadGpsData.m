% clc;clear all;close all;
function EphDat=ReadGpsData
global EphDat
EphDat=struct;
[filename,pathname]=uigetfile('*.**nav','��ȡGPS�㲥�����ļ�');%��nav�ļ���������filename��pathname���ļ�����·��
fid1=fopen(strcat(pathname,filename),'rt');%�򿪸��ļ�����ֻ��
%�ж��ļ���ȡ�Ƿ���ȷ
if(fid1==-1)
    msgbox('Input File or Path is not correct','warning ','warn');
   return;
end
%�ж���������
while(1)
    temp=fgets(fid1);
    header=findstr(temp,'END OF HEADER');%�ҡ�END OF HEADER�����ַ�λ�ò�����header
    if(~isempty(header))%�ж��Ƿ�Ϊ��
        break;
    end
end
i=1;
%��ʼд��ṹ���ļ�
while(1)
        temp=fgets(fid1)%ÿ��һ��������һ��
        if(temp==-1)
            break;
        end
        EphDat(i).PRN=str2num(temp(1:2));%��¼���ǵ�PRN��
        year=str2num(temp(4:5));%�����Ϣ
        EphDat(i).toc=TimetoJD(year,str2num(temp(7:8)),str2num(temp(10:11)),str2num(temp(13:14)),str2num(temp(16:17)),str2num(temp(19:22)));%����ʱ����
        EphDat(i).tocG=epoch2time(year,str2num(temp(7:8)),str2num(temp(10:11)),str2num(temp(13:14)),str2num(temp(16:17)),str2num(temp(19:22)));%ת��Ϊ����
        EphDat(i).a0=str2num(temp(23:41)); 
        EphDat(i).a1=str2num(temp(42:60));
        EphDat(i).a2=str2num(temp(61:79));
    temp=fgets(fid1);
        EphDat(i).idoe=str2num(temp(4:22));
        EphDat(i).Crs=str2num(temp(23:41));
        EphDat(i).dn=str2num(temp(42:60));
        EphDat(i).M0=str2num(temp(61:79));
    temp=fgets(fid1);
        EphDat(i).Cuc=str2num(temp(4:22));
        EphDat(i).e=str2num(temp(23:41));
        EphDat(i).Cus=str2num(temp(42:60));
        EphDat(i).sqa=str2num(temp(61:79));
    temp=fgets(fid1);
        EphDat(i).toes=str2num(temp(4:22));
        EphDat(i).Cic=str2num(temp(23:41));
        EphDat(i).omg0=str2num(temp(42:60));
        EphDat(i).Cis=str2num(temp(61:79));
    temp=fgets(fid1);
        EphDat(i).i0=str2num(temp(4:22));
        EphDat(i).Crc=str2num(temp(23:41));
        EphDat(i).w=str2num(temp(42:60));
        EphDat(i).omg=str2num(temp(61:79));
    temp=fgets(fid1);
        EphDat(i).iDot=str2num(temp(4:22));
        EphDat(i).cflgl2=str2num(temp(23:41));
        EphDat(i).weekno=str2num(temp(42:60));
        EphDat(i).pflgl2=str2num(temp(61:79));
        %EphDat(i).toe=adjweek(gpst2time(EphDat(i).weekno,EphDat(i).toes),EphDat(i).tocG);
        EphDat(i).toe=EphDat(i).tocG;
    temp=fgets(fid1);
        EphDat(i).svacc=uraindex(str2num(temp(4:22)));
        EphDat(i).svhlth=str2num(temp(23:41));
        EphDat(i).tgd=str2num(temp(42:60));
        EphDat(i).aodc=str2num(temp(61:79));
    temp=fgets(fid1);
        EphDat(i).ttm=str2num(temp(4:22));
    if(i~=1)                 %ɾ���ظ�����
         for k= i-1:i
             if (EphDat(i-1).PRN~=EphDat(i).PRN)
                 break;
             elseif abs(EphDat(i-1).toc-EphDat(i).toc)<0.001
                  i=i - 1;
             end
          end
     end
    i=i + 1;
end

function a=uraindex(value)

ura_eph=[ 2.4,3.4,4.85,6.85,9.65,13.65,24.0,48.0,96.0,192.0,384.0,768.0,1536.0,3072.0,6144.0,0.0];
for i=1:15
    if ura_eph(i)>=value
        a=i;
        break;
    end
end
    
