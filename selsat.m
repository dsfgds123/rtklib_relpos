function [ns,sat,ir,ib]=selsat(ObsODat_rover_k,nr,ObsODat_base_k,nb,azel_b)
% input
% ObsODat_rover_k :     移动站卫星
% nr :                  移动站卫星数
% ObsODat_base_k  :     基站卫星
% nb :                  基站卫星数
% azel_b ：             基站仰角方位角
% **************************
% output
% ns ：                 共视卫星数            
% sat ：                共视卫星号
% ir ：                 移动站共视卫星
% ib ：                 基站共视卫星
% *************************
ns=1;
i=1;
j=1;

for i=1:nr
    for j=1:nb
        if ObsODat_rover_k.SatCode(i)==ObsODat_base_k.SatCode(j)&azel_b(2*j)>=0.261799387799149 %基站仰角过低，信号不稳定，可视为该卫星质量不好
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
