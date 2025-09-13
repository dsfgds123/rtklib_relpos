function sol=pointpos(solrr,ObsODat_rover_k,EphDat,N)
% *************************************************************************
% input:
% solrr                        ��һ��Ԫ�Ľ��ջ�λ�ã���һ��ԪΪ0
% ObsODat_rover_k              ����Ԫ�Ľ��ջ��۲�ֵ
% EphDat                       ��������
% N                            ���������е����Ǹ����������ظ��ģ�
% *************************************************************************
% output:
% sol                          ���㶨λ�Ľ�
% *************************************************************************
%ֻ�ܶ�ȡU-blox��ʽ�ĵ�GPS����

tr=ObsODat_rover_k.TimeOEp;                %��Ԫ��������ʱ�䣨����Ԫ���൱��obs[0].time��
Snum=ObsODat_rover_k.SatSum;               %����Ԫ�۲⵽����������N��Sat.Sum��ͬ��
if Snum<4 
    return;                                %ȥ���޽����Ԫ
end
Code=ObsODat_rover_k.SatCode;              %����Ԫ�۲⵽�����Ǻ�
R=ObsODat_rover_k.Obs_RangeC1;             %ȡ��C1�۲�ֵ,��������ÿ�ֹ۲���������һ���У�

time=epoch2time(tr(1),tr(2),tr(3),tr(4),tr(5),tr(6));%�źŽ���ʱ��������ʱ����
tJD=TimetoJD(tr(1),tr(2),tr(3),tr(4),tr(5),tr(6));%������ʱ��
sol.time=time;%ʱ�����Ϊ�źŽ���ʱ��

%����λ��
[rs,dts,vare]=satposs(sol.time,tJD,ObsODat_rover_k,Snum,EphDat,R,N);

%���ջ�λ��
sol=estpos(time,solrr,rs,dts,Snum,R,Code,vare,EphDat);
