function [ xk1k,pk1k ] = rgt4( xkk ,pkk ,uk, Q, h)
    
    F = [-1 1; -0.2*xkk(1) 0];

    xp = [-xkk(1)+xkk(2); -0.1*xkk(1)^2 - 1 + uk];
    Pp = F*pkk + pkk*F' + Q;
    k1x = h*xp;
    k1p = h*Pp;
    
    xp = [-xkk(1)- k1x(1)/2 + xkk(2) + k1x(2)/2; -0.1*(xkk(1)+k1x(1)/2)^2 - 1 + uk];
    Pp = F*(pkk + k1p/2) + (pkk*k1p/2)*F' + Q;
    k2x = h*xp;
    k2p = h*Pp;
    
    xp = [-xkk(1)- k2x(1)/2 + xkk(2) + k2x(2)/2; -0.1*(xkk(1)+k2x(1)/2)^2 - 1 + uk];
    Pp = F*(pkk + k2p/2) + (pkk*k2p/2)*F' + Q;
    k3x = h*xp;
    k3p = h*Pp;
    
    xp = [-xkk(1)- k3x(1) + xkk(2) + k3x(2); -0.1*(xkk(1)+k3x(1))^2 - 1 + uk];
    Pp = F*(pkk + k3p) + (pkk*k3p)*F' + Q;
    k4x = h*xp;
    k4p = h*Pp;

    xk1k = xkk + k1x/6 + k2x/3 + k2x/3 + k4x/6;
    pk1k = pkk + k1p/6 + k2p/3 + k2p/3 + k4p/6;
    
end

