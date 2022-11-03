function [time,states] = simulate_quadrotor(fp,ts,x0,zeta0,xi0)

t0 = ts(1); tf = ts(2);

time   = zeros(length(t0:fp.T:tf),1);
states = zeros(length(t0:fp.T:tf),14);

t = t0;
X = [x0(1:9);zeta0;xi0;x0(10:12)];
% [x y z psi theta phi u v w zeta xi p q r]

time(1)     = t;
states(1,:) = X(:);

k = 1;
while t<=tf
    
    t = t+fp.T;

    Y = measurements(X(:),fp);
    
    [xref,yref,zref,psiref] = reference_trajectory(t,fp);
    ref_traj = {xref,yref,zref,psiref};
    
    [u1,u2,u3,u4] = controller(t,X(:),ref_traj,fp);
    control_input = {u1,u2,u3,u4};
    
    disturbances = external_disturbances();
    
    % Simulate with Forward Euler Integration
    X = state_propagation(t,X(:),control_input,disturbances,fp);
    
    k = k+1;
    time(k)     = t;
    states(k,:) = X(:);
    
end