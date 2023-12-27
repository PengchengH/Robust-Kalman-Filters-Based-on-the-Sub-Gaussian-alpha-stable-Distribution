function [xpos,P_pos] = kftncm(xt_,Pt_,z,s,Qn,Rn,F,H)

    % correct the measurement covariance with real scale value
    R = s*Rn;
   
    % Implement the Kalman Filter
    x_pre = F *xt_;
    P_pre = F *Pt_*F'+ Qn;
    K     = P_pre*H'*(R+H*P_pre*H')^-1;
    xpos  = x_pre +K*(z-H*x_pre);
    P_pos = (eye(size(P_pre,1))-K*H)*P_pre;   

end

