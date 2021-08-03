clc
clear all
close all

x_ = zeros(2,1);
P_ = eye(2);
Ts = 0.1;
tf = 100;
Nf = tf/Ts;
R = 0.01;
Q = 0.01*eye(2);
Hk1 = [1 0];
x0 = x_ + sqrt(P_)*randn(2,1);

sim('simulacao_1')

xkk(:,1) = [1;0];
pkk(:,:,1) = P_;

for k = 1:Nf

    yk1 = Dadosy.signals.values(k+1);
    uk = Dadosu.signals.values(k);
    
    [xk1k,Pk1k] = rgt4(xkk(:,k),pkk(:,:,k),uk, Q, Ts);
    
    
    Kk1 = Pk1k*Hk1'*inv(Hk1*Pk1k*Hk1'+R);
    rk1 = yk1 - xk1k(1);
    xkk(:,k+1) = xk1k + Kk1*rk1;
    pkk(:,:,k+1) = (eye(2)-Kk1*Hk1)*Pk1k;
    
end

x = Dadosx.signals.values;
xtil = x - xkk';

figure
hold;
plot(xkk(1,:),'r')
plot(x(:,1),'b')

figure
hold;
plot(xkk(2,:),'r')
plot(x(:,2),'b')

for k=1:Nf

        s1k(k) = sqrt(pkk(1,1,k));
        s2k(k) = sqrt(pkk(2,2,k));
end

figure; grid;
hold;
plot(xtil(:,1),'b');
plot(s1k,'r');
plot(-s1k,'r');

figure; grid;
hold;
plot(xtil(:,2),'b');
plot(s2k,'r');
plot(-s2k,'r');
