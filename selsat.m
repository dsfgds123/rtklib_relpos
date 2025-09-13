function [ns,sat,ir,ib]=selsat(ObsODat_rover_k,nr,ObsODat_base_k,nb,azel_b)
% input
% ObsODat_rover_k :     �ƶ�վ����
% nr :                  �ƶ�վ������
% ObsODat_base_k  :     ��վ����
% nb :                  ��վ������
% azel_b ��             ��վ���Ƿ�λ��
% **************************
% output
% ns ��                 ����������            
% sat ��                �������Ǻ�
% ir ��                 �ƶ�վ��������
% ib ��                 ��վ��������
% *************************
ns=1;
i=1;
j=1;

for i=1:nr
    for j=1:nb
        if ObsODat_rover_k.SatCode(i)==ObsODat_base_k.SatCode(j)&azel_b(2*j)>=0.261799387799149 %��վ���ǹ��ͣ��źŲ��ȶ�������Ϊ��������������
            sat(ns)=ObsODat_rover_k.SatCode(i);
            ir(ns)=i;
            ib(ns)=j;
            ns=ns+1;
            break;
        end
    end
end
% while i<=nr&j<=nb
%     if ObsODat_rover_k.SatCode(i)<ObsODat_base_k.SatCode(j)
%         j=j-1;
%     elseif ObsODat_rover_k.SatCode(i)>ObsODat_base_k.SatCode(j)
%         i=i-1;
%     elseif azel_b(2*j)>=0.261799387799149
%         sat(ns)=ObsODat_rover_k.SatCode(i);
%         ir(ns)=i;
%         ib(ns)=j;
%         ns=ns+1;
%     end
%     i=i+1;
%     j=j+1;
% end
ns=ns-1;
