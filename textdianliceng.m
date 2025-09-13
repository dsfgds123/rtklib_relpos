t=[1653864349 0];
pos=[25.107 110.371];
azel=[292.376 33.728];
CLIGHT=299792458;
        ion_default=[0.1118E-07 -0.7451E-08 -0.5961E-07 0.1192E-06 0.1167E+06 -0.2294E+06 -0.1311E+06 0.1049E+07];
        if(azel(2)<=0)
            dion=0;
        else
            ion=ion_default;
            psi=0.0137/(azel(2)/pi+0.11)-0.022;
            phi=pos(1)/pi+psi*cos(azel(1));
            if(phi>0.416)
                phi=0.416;
            elseif(phi<-0.416)
                phi=-0.416;
            end
            LAM=pos(2)/pi+psi*sin(azel(1))/cos(phi*pi);
            phi=phi+0.064*cos((LAM-1.617)*pi);
            week=time2gpst(t);%用到输入量，信号接收时间
            tt=43200.0*LAM+week(2);
            tt=tt-floor(tt/86400.0)*86400.0;%0<=tt<86400.0
            f=1.0+16.0*(0.53-azel(2)/pi)^3;
            amp=ion(1)+phi*(ion(2)+phi*(ion(3)+phi*ion(4)));
            per=ion(5)+phi*(ion(6)+phi*(ion(7)+phi*ion(8)));
            if(amp<0)
                amp=0;
            end
            if(per<72000.0)
                per=72000.0;
            end
            x=2.0*pi*(tt-50400)/per;
            if(abs(x)<1.57)
                temp=5E-9+amp*(1.0+x*x*(-0.5+x*x/24.0));
            else
                temp=5E-9;
            end
            dion=CLIGHT*f*temp;%电离层修正
        end
        vion=(dion*0.5)*(dion*0.5);%电离层方差