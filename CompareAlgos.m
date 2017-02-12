clear all

inpath = 'report';
%inpath  = 'data/MATLAB-data';

G_orig = [
    0.10,0.56,0.30,0.89,0.02;
    0.50,0.10,0.15,0.05,0.20;
    0.25,0.30,0.10,0.10,0.25]';
G_orig = G_orig*diag(1./sum(G_orig,1));

N = 5;
K = 3;
Gs_data = csvread(strcat(inpath,'/PL_Gs_K',num2str(K),'_N',num2str(N),'-inferred-noprior.csv'));
[MSE_1, nepochs_1] = computeMSE(Gs_data, G_orig, N, K);

Gs_data = csvread(strcat(inpath,'/PL_Gs_K',num2str(K),'_N',num2str(N),'-inferred-intpoint.csv'));
[MSE_2, nepochs_2] = computeMSE(Gs_data, G_orig, N, K);

figure();
plot(1:nepochs_1, MSE_1, '--b');
hold on;
plot(1:nepochs_2, MSE_2, '-r');
hold off;
