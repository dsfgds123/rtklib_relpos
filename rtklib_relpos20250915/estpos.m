function [rr_new, dtr_new, stat_new]=estpos(t,solkrr,rs,dts,Snum,R,Code,vare,nav,Sys)
% *************************************************************************
% --- FINAL VERSION: Returns raw values instead of an incomplete struct ---
% *************************************************************************

global CLIGHT;

% --- Initialize outputs to a safe 'no solution' state ---
rr_new = [0,0,0];
dtr_new = 0;
stat_new = 0; % SOLQ_NONE

X=[solkrr(1);solkrr(2);solkrr(3);0];

 for i=1:10
     [v,var,nv,H]=PNTrescode(X,t,rs,dts,vare,Snum,R,Code,nav,Sys);
     
     if(nv<4), break; end % Not enough valid satellites
     
     for j=1:nv
         sig=sqrt(var(j));
         v(j)=v(j)/sig;
         H(j,:)=H(j,:)/sig;
     end
     
     if rcond(H'*H) < 1e-12, break; end % Avoid singular matrix
     dx=(H'*H)^-1*H'*v'; % Note: v is a column, H'*v needs v'
     
     X=X+dx;
     
     if(norm(dx)<1E-4)
         % --- MODIFIED: Return raw values ---
         rr_new = [X(1), X(2), X(3)];
         dtr_new = X(4);
         stat_new = 1; % SOLQ_SINGLE
         return; % Solution converged, exit function
     end
 end
end