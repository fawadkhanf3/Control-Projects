clear; clc;
FLAG = 1; % 0 for Ours, 1 for Boyd 

mdl = 'Ex01_eval';

load_system(mdl);

io(1) = linio('Ex01_eval/w',1,'input');
io(2) = linio('Ex01_eval/z',1,'output');

A = 1; %[0.01, 0.1, 1, 10];
w = logspace(-3,3,25);
tic
for ct = 1:numel(A)
    in = frest.Sinestream('Frequency',w,'Amplitude',A(ct),...
        'NumPeriods',10,'SettlingPeriods',7);
    G(1,1,ct)  = frestimate(mdl,in,io);
    disp(ct)
end
toc