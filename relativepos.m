function rtk=relativepos(rtk,ObsODat_rover_k,ObsODat_base_k,EphDat,N,mode,ARmode)

global MINFIX;
global SOLQ_NONE;%�޽�
global SOLQ_FLOAT;%�����
global SOLQ_FIX;%������
global MAXSAT;

stat=SOLQ_FLOAT;

Snum_r=ObsODat_rover_k.SatSum; Snum_b=ObsODat_base_k.SatSum;                 %��ǰ��Ԫ��վ�ƶ�վ��������

R_r=ObsODat_rover_k.Obs_RangeC1; R_b=ObsODat_base_k.Obs_RangeC1;             %��վ�ƶ�վ�����ǵ�α�����ֵ

tr_r=ObsODat_rover_k.TimeOEp;                                                  %��վ�ƶ�վ����Ԫʱ��
time_r=epoch2time(tr_r(1),tr_r(2),tr_r(3),tr_r(4),tr_r(5),tr_r(6));
tr_b=ObsODat_base_k.TimeOEp;
time_b=epoch2time(tr_b(1),tr_b(2),tr_b(3),tr_b(4),tr_b(5),tr_b(6));

dt=timediff(time_r,time_b);                                                    %��վ�ƶ�վ������������ʱ���

tJD_r=TimetoJD(tr_r(1),tr_r(2),tr_r(3),tr_r(4),tr_r(5),tr_r(6));               %������ʱ��

[RS_r,DTS_r,VAR_r]=satposs(time_r,tJD_r,ObsODat_rover_k,Snum_r,EphDat,R_r,N); %�ƶ�վ����λ��  
[RS_b,DTS_b,VAR_b]=satposs(time_r,tJD_r,ObsODat_base_k,Snum_b,EphDat,R_b,N);  %��վ����λ��-
 %��վ�ǲ�ֲв�
[y_b,e_b,azel_b]=zdres(time_b,ObsODat_base_k,Snum_b,RS_b,DTS_b,rtk.rb); %��һ�εķǲ�ֲв���Ϊ�����Ƿ�λ��ȥѰ�ҹ�������                 

%time-interpolation of residuals (for post-processing)
%dt=intpres(time_r,ObsODat_base_k,Snum_b,EphDat,rtk,y_b);

[ns,sat,ir,ib]=selsat(ObsODat_rover_k,Snum_r,ObsODat_base_k,Snum_b,azel_b);%���ҹ�������

rtk=udstate(mode,ARmode,rtk,ObsODat_rover_k,ObsODat_base_k,sat,ir,ib,ns,EphDat);%����RTK�ṹ�壨��ǰ��Ԫ��ʼ����

% if abs(2005009.6688215812+rtk.x(1))>0.1 %��һ���жϣ�Ӧ���Ǻ���ʵ������һ���Ƚ�
%     aaa=1;
%     bbb=2;   %�������� 
% end

xp=rtk.x;%���ջ�״̬
% 1�ǵ����������ɸ�����Ҫ������������
% for i=1:1                                    %���ѭ����ɶ��˼����
    %�ƶ�վ�ǲ�ֲв�
    [y_r,e_r,azel_r]=zdres(time_r,ObsODat_rover_k,Snum_r,RS_r,DTS_r,xp);%����ز���λα��в�ǵص�λʸ���Լ����Ƿ�λ��
    
    for j=1:Snum_r
        for k=1:MAXSAT
            if ObsODat_rover_k.SatCode(j)==k
                rtk.ssat(k).azel=[azel_r(2*j-1),azel_r(2*j)];%��������λ��  �ѽ�������ķ�λ���������ݴ���RTK��
                break;
            end
        end
    end

    %���Ĳ���1����ddres�и�H��ֵ��0������ֵ
    [nv,rtk,v,H,R,vflg]=ddres(rtk,dt,xp,sat,y_r,e_r,azel_r,y_b,ir,ib,ns,1);%%˫�����
    
    Pp=rtk.P;%%����Э������
    
    [xp_,Pp_]=EKF_filter(xp,Pp,H,v,R,rtk.nx,nv,0);%�������˲����˳����λ������ջ�����
% end

if(stat~=SOLQ_NONE) %�������ⲻ���޽�                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                  
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                    
    [y_r,e_r,azel_r]=zdres(time_r,ObsODat_rover_k,Snum_r,RS_r,DTS_r,xp_);%�������˲������ٽ��зǲ����
    
    %�˴�H��ֵ����EKF_filter�õ�ֵ��ͬ��û�б仯
    [nv,rtk,v,H,R,vflg]=ddres(rtk,dt,xp_,sat,y_r,e_r,azel_r,y_b,ir,ib,ns,0);%�˲�֮��Ĳв��˫��
    
    if(valpos(rtk,v,R,vflg,nv,4))%��Ч�Լ���
        %����״̬��Э�������
        rtk.x=xp_;
        rtk.P=Pp_;
        %����ģ���ȿ��ƽṹ��
        rtk.sol.ns=0;%��Ч������
        for i=1:ns             %%����������
            if ~rtk.ssat(sat(i)).vsat%%�����Ƿ���Ч��־
                continue;
            end
            rtk.ssat(sat(i)).lock=rtk.ssat(sat(i)).lock+1;%��λ��������
            rtk.ssat(sat(i)).outc=0;%�����źŶϹ�����
            rtk.sol.ns=rtk.sol.ns+1;
        end
        %ȱ����Ч����
        if rtk.sol.ns<4
            stat=SOLQ_NONE;
        end
    else
        stat=SOLQ_NONE;
    end
end

if stat~=SOLQ_NONE
    % [rtk,xa,AmbNum]=resamb_LAMBDA(rtk);
    [rtk,xa,AmbNum]=resamb_AGGWO(rtk);%�������ģ���ȼ���Ӧ���ջ�״̬
    if AmbNum>1
        [y_r,e_r,azel_r]=zdres(time_r,ObsODat_rover_k,Snum_r,RS_r,DTS_r,xa);%���º��xa����в�
        
        %����Ҳ�������H��������ddres�����Ĵ��ݲ��������һ������H_YorN����������H
        %�˴�Ҳ�����ô���Pp���������һ������Pp_YorN
        [nv,rtk,v,H,R,vflg]=ddres(rtk,dt,xa,sat,y_r,e_r,azel_r,y_b,ir,ib,ns,0);%���º����ֵ����˫��
        
        if(valpos(rtk,v,R,vflg,nv,4))%��Ч�Լ���
            %��������ģ����
            rtk.nfix=rtk.nfix+1;%�̶��������1
            if(rtk.nfix>=MINFIX&ARmode==2)
                [SsatFix,rtk.x,rtk.P]=holdamb(rtk,xa);
                for i=1:MAXSAT
                    rtk.ssat(i).fix=SsatFix(i);%����״̬�ṹ��̶���
                end
            end
            stat=SOLQ_FIX;
        end
    end
end
%�����
if stat==SOLQ_FIX
    rtk.sol.rr=[rtk.xa(1), rtk.xa(2), rtk.xa(3)];
    rtk.sol.qr=[rtk.Pa(1,1), rtk.Pa(2,2), rtk.Pa(3,3)];
else
    rtk.sol.rr=[rtk.x(1), rtk.x(2), rtk.x(3)];
    rtk.sol.qr=[rtk.P(1,1), rtk.P(2,2), rtk.P(3,3)];
    rtk.nfix=0;
end

for i=1:MAXSAT
    if rtk.ssat(i).fix==2&stat~=SOLQ_FIX       %�������״̬�̶���Ϊ2���Ҹ���ⲻ�̶�
        rtk.ssat(i).fix=1;                     %�����Ϊ1��
    end
end
