function [head,obs]=ReadObsData
%读接收机观测数据文件*******************************************************
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

        [fname,fpath]=uigetfile('*.**O','选择一个O文件');
        O_filename=strcat(fpath,fname);
   
    fid1=fopen(O_filename,'rt');
    if (fid1==-1)
        msgbox('file invalide','warning','warn');
        return;
    end
    %将文件头数据存入结构体HeadODat中
    t=0;
    while(t<100)
        s=fgets(fid1);
        t=t+1;
        L=size(s,2);
        if L<81
            s(L+1:81)=' ';
        end
          
        instrS=s(61:81);
        %测站点近似坐标
        if strncmp(instrS,'APPROX POSITION XYZ',19)
            HeadODat.ApproXYZ=zeros(1,3);
            HeadODat.ApproXYZ(1,1)=str2num(s(1:14));
            HeadODat.ApproXYZ(1,2)=str2num(s(15:28));
            HeadODat.ApproXYZ(1,3)=str2num(s(29:42));
        %历元间隔;
        elseif strncmp(instrS,'INTERVAL',8)
            HeadODat.interval=str2num(s(5:11));
        %测站点号    
        elseif strncmp(instrS,'MARKER NAME',11)
            HeadODat.SiteName=s(1:4)
        %天线中心改正    
        elseif strncmp(instrS,'ANTENNA: DELTA H/E/N',20)
            HeadODat.Ant_H=str2num(s(1:14));
            HeadODat.Ant_E=str2num(s(15:28));
            HeadODat.Ant_N=str2num(s(29:42));
        %第一个历元时间    
        elseif strncmp(instrS,'TIME OF FIRST OBS',17)
            year=str2num(s(1:6));
            month=str2num(s(7:12));
            day=str2num(s(13:18));
            hour=str2num(s(19:24));
            minute=str2num(s(25:30));
            second=str2num(s(31:42));
            HeadODat.TimeOB=TimetoJD(year-2000,month,day,hour,minute,second);
        %最后一个历元时间    
        elseif strncmp(instrS,'TIME OF LAST OBS',16)
            year=str2num(s(1:6));
            month=str2num(s(7:12));
            day=str2num(s(13:18));
            hour=str2num(s(19:24));
            minute=str2num(s(25:30));
            second=str2num(s(31:42));
            HeadODat.TimeOE=TimetoJD(year-2000,month,day,hour,minute,second);
        %观测值类型    
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
        %头文件结束  
        
        
        elseif strncmp(instrS,'END OF HEADER',13)
            break;
       else
            continue;
       end 
   end

   
   
   %观测数据结构体%观测数据结构
        t=0;
   while ~feof(fid1)
        %每个历元的第一行数据，时间和观测到的卫星号
        s=fgets(fid1);
        t=t+1;
        year=str2num(s(1:3));
        month=str2num(s(4:6));
        day=str2num(s(7:9));
        hour=str2num(s(10:12));
        minute=str2num(s(13:15));
        second=str2num(s(16:26));
        %历元时间 
        ObsODat(t).TimeOEp=[year,month,day,hour,minute,second];
        
        ObsODat(t).TimeOEpp=TimetoJD(year,month,day,hour,minute,second);
        %该历元观测到的卫星数
        ObsODat(t).SatSum=str2num(s(30:32));
        %该历元观测到的卫星号
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
       
        %每个历元的观测数据，按卫星号先后顺序分行存
        for k=1:ObsODat(t).SatSum
            s=fgets(fid1);
            %判断一个卫星的观测数据是否分两行存储，如果为两行，则再读入一行
            if HeadODat.SumOType>5
                sg=fgets(fid1);
                s=strcat(s,sg);
            end
            L=size(s,2);
            %补充数据长度
            if L<HeadODat.SumOType*16
                s(L+1:HeadODat.SumOType*16)=' ';
            end

        
            %对观测数据判断其类型，并存储到相应的数组中
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
            %完成一个卫星的所有观测数据存储
        end
        %完成一个历元的数据存储
    end 
    %完成所有历元的数据存储
    
    head=HeadODat;
    obs=ObsODat;
    return 