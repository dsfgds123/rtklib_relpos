function [head,obs]=ReadObsData(a)
% version history *********************************************************
% 2017/05/12
% ������һ��if���"if a==0"������ʾ˵�������ƶ�վ���վ����
% ������һ�����"HeadODat.TimeStart=[year,month,day,hour,minute,second];"
% ������һ�����"HeadODat.TimeEnd=[year,month,day,hour,minute,second];"
% 2017/05/17
% �������ź�ʧ����־LLI��ֻ������ز��۲���L1���ź�ʧ����־
% ���Ϊ��ObsODat(t).LLI(1,k)=stempNumLLI;
% 2017/05/24
% ���������飬ʹ������Զ��ǿ�������12����Щ��Ԫ
% ���Ϊ if(ObsODat(t).SatSum <= 12)�� if(ObsODat(t).SatSum > 12)
% version history *********************************************************


%�����ջ��۲������ļ�*******************************************************
% HeadODat :a structure stores header information if o-file
% .ApproXYZ[3];              approximate coordinate����վ�������꣩
% .interval;                 intervals of two adjacent epochs�����������ڵļ����
% .SiteName;                 point name��վ�����ƣ�
% .Ant_H;                    antenna height�����߸߶ȣ�
% .Ant_E;                    east offset of the antenna center�����߶���ƫ�ģ�
% .Ant_N;                    north offset of then antenna center�����߱���ƫ�ģ�
% .TimeOB;                   first epoch time in modifuied Julian time����һ�������ļ�ʱ�䣬����������ʱ���б�ʾ��
% .TimeOE;                   last epoch time in modifuied Julian time�����һ�������ļ�ʱ�䣬����������ʱ���б�ʾ��
% .SumOType;                 number of types of observables���ɹ۲�������͵�������
% .SumOO[SumOType];          type of observables
% 0-L1,1-L2,2-C1,3-P1,4-P2,5-D1,6-D2,���ɹ۲�������ͣ�

%ObsODat :a structure stores observables by one and one epoch
% .TimeOEpp[2];              reciever time of epoch
% .SatSum;                   number of satellites
% .SatCode[SatSum];          satellites' PRN
% .Obs_FreL1[SatSum]; 
% .Obs_FreL2[SatSum];
% .Obs_RangeC1[SatSum];
% .Obs_RangeP1[SatSum];
% .Obs_RangeP2[SatSum];
% *************************************************************************
    if a==0
        [fname,fpath]=uigetfile('*.*','ѡ���ƶ�վO�ļ�');
    else
        [fname,fpath]=uigetfile('*.*','ѡ���վO�ļ�');
    end
    O_filename=strcat(fpath,fname);
   
    fid1=fopen(O_filename,'rt');
    if (fid1==-1)
        msgbox('file invalide','warning','warn');
        return;
    end
    %���ļ�ͷ���ݴ���ṹ��HeadODat��
    t=0;
    while(t<100)
        s=fgets(fid1);
        t=t+1;%t��־ͷ�ļ�������
        L=size(s,2);
        if L<81
            s(L+1:81)=' ';%�ÿո����һ����û�����ݵĵط���ʹÿһ�ж���81���ַ�
        end
          
        instrS=s(61:81);
        %��վ���������
        if strncmp(instrS,'APPROX POSITION XYZ',19)
            HeadODat.ApproXYZ=zeros(1,3);
            HeadODat.ApproXYZ(1,1)=str2num(s(1:14));
            HeadODat.ApproXYZ(1,2)=str2num(s(15:28));
            HeadODat.ApproXYZ(1,3)=str2num(s(29:42));
        %��Ԫ���;
        elseif strncmp(instrS,'INTERVAL',8)
            HeadODat.interval=str2num(s(5:11));
        %��վ���
        elseif strncmp(instrS,'MARKER NAME',11)
            HeadODat.SiteName=s(1:4)
        %�������ĸ���
        elseif strncmp(instrS,'ANTENNA: DELTA H/E/N',20)
            HeadODat.Ant_H=str2num(s(1:14));
            HeadODat.Ant_E=str2num(s(15:28));
            HeadODat.Ant_N=str2num(s(29:42));
        %��һ����Ԫʱ��    
        elseif strncmp(instrS,'TIME OF FIRST OBS',17)
            year=str2num(s(1:6));
            month=str2num(s(7:12));
            day=str2num(s(13:18));
            hour=str2num(s(19:24));
            minute=str2num(s(25:30));
            second=str2num(s(31:42));
            HeadODat.TimeStart=[year,month,day,hour,minute,second];
            HeadODat.TimeOB=TimetoJD(year-2000,month,day,hour,minute,second);
        %���һ����Ԫʱ��    
        elseif strncmp(instrS,'TIME OF LAST OBS',16)
            year=str2num(s(1:6));
            month=str2num(s(7:12));
            day=str2num(s(13:18));
            hour=str2num(s(19:24));
            minute=str2num(s(25:30));
            second=str2num(s(31:42));
            HeadODat.TimeEnd=[year,month,day,hour,minute,second];
            HeadODat.TimeOE=TimetoJD(year-2000,month,day,hour,minute,second);
        %�۲�ֵ����    
        elseif strncmp(instrS,'# / TYPES OF OBSERV',19)
            HeadODat.SumOType=str2num(s(1:6));%�۲��������ͣ�α�ࡢ�ز���λ����S1
            HeadODat.SumOO=ones(1,HeadODat.SumOType)*-1;
            for k=1:HeadODat.SumOType
                f=s(k*6+5:k*6+6); 
               x=strcmp(f,'C1')
                if strcmp(f,'L1')
                    HeadODat.SumOO(1,k)=0;
                elseif strcmp(f,'L2')
                    HeadODat.SumOO(1,k)=1;
                elseif strcmp(f,'C1')
                        HeadODat.SumOO(1,k)=2;
                elseif strcmp(f,'P1')
                        HeadODat.SumOO(1,k)=3;
                elseif strcmp(f,'P2')
                        HeadODat.SumOO(1,k)=4;
                elseif strcmp(f,'D1')
                        HeadODat.SumOO(1,k)=5;
                elseif strcmp(f,'D2')
                        HeadODat.SumOO(1,k)=6;
                end
            end
        %ͷ�ļ�����  
       elseif strncmp(instrS,'END OF HEADER',13)
            break;
       else
            continue;
       end 
   end

   
   
   %�۲����ݽṹ��%�۲����ݽṹ
        t=0;
   while ~feof(fid1)
        %ÿ����Ԫ�ĵ�һ�����ݣ�ʱ��͹۲⵽�����Ǻ�
        s=fgets(fid1);
        t=t+1;
        year=str2num(s(1:3));
        month=str2num(s(4:6));
        day=str2num(s(7:9));
        hour=str2num(s(10:12));
        minute=str2num(s(13:15));
        second=str2num(s(16:26));
        %��Ԫʱ�� 
        ObsODat(t).TimeOEp=[year,month,day,hour,minute,second];
        
        ObsODat(t).TimeOEpp=TimetoJD(year,month,day,hour,minute,second);
        %����Ԫ�۲⵽��������
        ObsODat(t).SatSum=str2num(s(30:32));
        %����Ԫ�۲⵽�����Ǻ�
        ObsODat(t).SatCode=zeros(1,ObsODat(t).SatSum);
        ObsODat(t).Obs_FreL1=zeros(1,ObsODat(t).SatSum);
        ObsODat(t).Obs_FreL2=zeros(1,ObsODat(t).SatSum);
        ObsODat(t).Obs_RangeC1=zeros(1,ObsODat(t).SatSum);
        ObsODat(t).Obs_RangeP1=zeros(1,ObsODat(t).SatSum);
        ObsODat(t).Obs_RangeP2=zeros(1,ObsODat(t).SatSum);
        ObsODat(t).LLI=zeros(1,ObsODat(t).SatSum);
        if(ObsODat(t).SatSum <= 12)
            for k=1:ObsODat(t).SatSum
                f=s(31+k*3:32+k*3);
                ObsODat(t).SatCode(1,k)=str2num(f);
            end        
        end
        if (ObsODat(t).SatSum > 12)
            for k=1:12
                f=s(31+k*3:32+k*3);
                ObsODat(t).SatCode(1,k)=str2num(f);
            end
            s=fgets(fid1);
            for k=1:ObsODat(t).SatSum-12
                f=s(31+k*3:32+k*3);
                ObsODat(t).SatCode(1,k+12)=str2num(f);
            end            
        end
       
        %ÿ����Ԫ�Ĺ۲����ݣ������Ǻ��Ⱥ�˳����д�
        for k=1:ObsODat(t).SatSum
            s=fgets(fid1);
            %�ж�һ�����ǵĹ۲������Ƿ�����д洢�����Ϊ���У����ٶ���һ��
            if HeadODat.SumOType>5
                sg=fgets(fid1);
                s=strcat(s,sg);
            end
            L=size(s,2);%����s���е��������ַ���
            %�������ݳ���
            if L<HeadODat.SumOType*16
                s(L+1:HeadODat.SumOType*16)=' ';%�Ѹ����ǺŶ�Ӧ��һ����û���ݵĵط��ÿո����
            end

        
            %�Թ۲������ж������ͣ����洢����Ӧ��������
            for j=1:HeadODat.SumOType
                stemp=s(j*16-15:j*16-2);% stemp �õ�����һ�����͵�һ�����ݱ���L1�ģ�     -6794.564 7
                stempLLI=s(j*16);%�õ��ź�ʧ����־
                stemp=deblank(stemp);%ֻɾ��β���Ŀո�
                stempLLI=deblank(stempLLI);
                if isempty(stemp)
                    continue;
                end
                if isempty(stempLLI)
                    stempNumLLI=0;
                else
                    stempNumLLI=str2num(stempLLI);
                end
                stempNum=str2num(stemp);
                stempLength=size(stempNum,2);
                if stempLength>1
                    stempNum=stempNum(1,1);
                end
                if HeadODat.SumOO(1,j)==0
                    ObsODat(t).Obs_FreL1(1,k)=stempNum;
                    ObsODat(t).LLI(1,k)=stempNumLLI;%ֻ�����λ��ʧ����־
                elseif HeadODat.SumOO(1,j)==1
                    ObsODat(t).Obs_FreL2(1,k)=stempNum;
                elseif HeadODat.SumOO(1,j)==2
                    ObsODat(t).Obs_RangeC1(1,k)=stempNum;
                elseif HeadODat.SumOO(1,j)==3
                    ObsODat(t).Obs_RangeP1(1,k)=stempNum;
                elseif HeadODat.SumOO(1,j)==4
                    ObsODat(t).Obs_RangeP2(1,k)=stempNum;
                else
                    continue;
                end
            end
            %���һ�����ǵ����й۲����ݴ洢
        end
        %���һ����Ԫ�����ݴ洢
   end 
    %���������Ԫ�����ݴ洢
    
    head=HeadODat;
    obs=ObsODat;
    return 