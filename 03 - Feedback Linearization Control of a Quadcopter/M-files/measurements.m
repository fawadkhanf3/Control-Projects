function Y = measurements(X,fp)

position_noise = fp.position_noise*(rand(2,1)*2-1);
height_noise   = fp.height_noise*(rand(1,1)*2-1);
angle_noise    = fp.angle_noise*(rand(3,1)*2-1);
rate_noise     = fp.rate_noise*(rand(3,1)*2-1);

xhat = X(1) + position_noise(1,1);
yhat = X(2) + position_noise(2,1);
zhat = X(3) + height_noise(1,1);

position_hat = [xhat;yhat;zhat];

psihat   = X(7) + angle_noise(1,1);
thetahat = X(8) + angle_noise(2,1);
phihat   = X(9) + angle_noise(3,1);

angle_hat = [psihat;thetahat;phihat];

phat = X(12) + rate_noise(1,1);
qhat = X(13) + rate_noise(2,1);
rhat = X(14) + rate_noise(3,1);

rate_hat = [phat;qhat;rhat];

Y = [position_hat;X(4:6);angle_hat;X(10:11);rate_hat];

end