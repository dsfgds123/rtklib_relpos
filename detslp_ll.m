function SLIP=detslp_ll(rtk,ObsODat_k,i,base)

sat=ObsODat_k.SatCode(i);

LLI1=bitand(bitshift(rtk.ssat(sat).slip,-6),3);
LLI2=bitand(bitshift(rtk.ssat(sat).slip,-4),3);

if base==1
    LLI=LLI1;
elseif base==2
    LLI=LLI2;
end

slip=bitand(bitor(rtk.ssat(sat).slip,ObsODat_k.LLI(i)),3);

if (((bitand(LLI,2))&~(bitand(ObsODat_k.LLI(i),2)))|(~(bitand(LLI,2))&(bitand(ObsODat_k.LLI(i),2))))
    slip=bitor(slip,1);
end

if base==1
    SLIP=bitor(bitshift(ObsODat_k.LLI(i),6),bitshift(LLI2,4));
    SLIP=bitor(SLIP,slip);
else
    SLIP=bitor(bitshift(ObsODat_k.LLI(i),4),bitshift(LLI1,6));
    SLIP=bitor(SLIP,slip);
end