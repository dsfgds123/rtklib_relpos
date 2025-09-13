function stat=valpos(rtk,v,R,vflg,nv,thres)
% *************************************************************************
% input:
% R                            ˫��������ķ���Э�������
% vflg                         ��־λ
% nv                           ˫��в����
% v                            ˫��в�����
% rtk                          ��������ṹ��
% *************************************************************************
% output:
% stat                         ?
% *************************************************************************
fact=thres*thres;
stat=1;

for i=1:nv
    if v(i)*v(i)<=fact*R(i,i)
        continue;
    end
    
    sat1=bitand(bitshift(vflg(i),-16),255);
    sat2=bitand(bitshift(vflg(i), -8),255);
    type=bitand(bitshift(vflg(i), -4), 15);
    freq=bitand(vflg(i),15);
    if type==0|type==1
        stype='L';
    else
        stype='C';
    end
end
    