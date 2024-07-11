%% 
clear all 
clc
load times_090425blk10_ch66.mat
%%
% Example spike trains for two neurons
T1=find(cluster_class(:,1)==1);
T1=T1/par.sr;
T2=find(cluster_class(:,1)==2);
T2=T2/par.sr;
max_length = max(length(T1), length(T2));
T2 = interp1(1:length(T2), T2, linspace(1, length(T2), max_length));
lags = -(length(T1)-1):(length(T2) -1);
cross_correlation = xcorr(T1,T2);
stem(lags, cross_correlation);
xlabel('Lag');
ylabel('Cross-correlation');
title('Cross-correlation between Neuron 1 and Neuron 2');

% Find the lag with maximum cross-correlation
[max_correlation, max_index] = max(cross_correlation);
lag_at_max_correlation = lags(max_index);

disp(['Maximum cross-correlation: ', num2str(max_correlation)]);
disp(['Lag at maximum cross-correlation: ', num2str(lag_at_max_correlation)]);
