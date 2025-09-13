function [head,obs]=ReadObsData
%�����ջ��۲������ļ�*******************************************************
% HeadODat :a structure stores header information if o-file
% .ApproXYZ[3];              approximate coordinate
% .interval;                 intervals of two adjacent epochs
% .SiteName;                 point name
% .Ant_H;                    antenna height
% .Ant_E;                    east offset of the antenna center
% .Ant_N;                    north offset of then antenna center
% .TimeOB;                   first epoch time in modifuied Julian time
% .TimeOE;                   last epoch time in modifuied Julian time
% .SumOType;                 number of types of observables
% .SumOO[SumOType];          type of observables 0-L1,1-L2,2-C1,3-P1,4-P2,5-D1,6-D2,

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
global HeadODat;
global ObsODat;

        [fname,fpath]=uigetfile('*.**O','ѡ��һ��O�ļ�');
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
        t=t+1;
        L=size(s,2);
        if L<81
            s(L+1:81)=' ';
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
            HeadODat.TimeOB=TimetoJD(year-2000,month,day,hour,minute,second);
        %���һ����Ԫʱ��    
        elseif strncmp(instrS,'TIME OF LAST OBS',16)
            year=str2num(s(1:6));
            month=str2num(s(7:12));
            day=str2num(s(13:18));
            hour=str2num(s(19:24));
            minute=str2num(s(25:30));
            second=str2num(s(31:42));
            HeadODat.TimeOE=TimetoJD(year-2000,month,day,hour,minute,second);
        %�۲�ֵ����    
        elseif strncmp(instrS,'# / TYPES OF OBSERV',19)
            HeadODat.SumOType=str2num(s(1:6));
            HeadODat.SumOO=ones(1,HeadODat.SumOType)*-1;
            for k=1:HeadODat.SumOType
                f=s(k*6+5:k*6+6); 
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
            for k=13:ObsODat(t).SatSum-12
                f=s(31+k*3:32+k*3);
                ObsODat(t).SatCode(1,k)=str2num(f);
            end            
        end

        %ObsODat(t).TDiff = str2num(s(70:81));
       
        %ÿ����Ԫ�Ĺ۲����ݣ������Ǻ��Ⱥ�˳����д�
        for k=1:ObsODat(t).SatSum
            s=fgets(fid1);
            %�ж�һ�����ǵĹ۲������Ƿ�����д洢�����Ϊ���У����ٶ���һ��
            if HeadODat.SumOType>5
                sg=fgets(fid1);
                s=strcat(s,sg);
            end
            L=size(s,2);
            %�������ݳ���
            if L<HeadODat.SumOType*16
                s(L+1:HeadODat.SumOType*16)=' ';
            end

        
            %�Թ۲������ж������ͣ����洢����Ӧ��������
            for j=1:HeadODat.SumOType
                stemp=s(j*16-15:j*16);
                stemp=deblank(stemp);
                if isempty(stemp)
              
                    
                    continue;
                end
                stempNum=str2num(stemp);
                stempLength=size(stempNum,2);
                if stempLength>1 
                    stempNum=stempNum(1,1);
                end
                if HeadODat.SumOO(1,j)==0
                    ObsODat(t).Obs_FreL1(1,k)=stempNum;
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