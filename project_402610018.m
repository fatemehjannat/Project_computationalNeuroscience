%% Fatemeh Jannat Fereidouni(402610018)
clear all
close all
clc
%% load data
load 'times_090425blk10_ch66.mat'
%% single trial analysis
par.scales=4;
data=spikes(:,12)*par.scales;
figure;
plot(data)
title('Raw data');
%% fit sampling rate in PSD
pwelch(data,[],[],[],par.sr);
%% filtering
order=3;
cutoff_freq=1100;%Hz
[b,a] = butter(order,cutoff_freq/(par.sr/2),'high');
filtered_data=filtfilt(b,a,data);
figure;
plot(filtered_data)
title('filtered data')
hold on
% spike detection
sig=median(abs(filtered_data)/0.6745);
plot(ones(length(filtered_data),1)*5*sig,'LineWidth',2)
spike_train=filtered_data>5*sig;
data_tmp=filtered_data;
data_tmp(filtered_data<5*sig)=0;
[spike,spike_index]=findpeaks(data_tmp,'MinPeakDistance',100);
%%
maxtime=length(data)/par.sr;
plot(linspace(0,maxtime,length(data)),data)
title('Neural Response Recording')
xlabel('time(s)')
ylabel('amplitude(mv)')
%% created spike train
spikeTrain=zeros(1,length(data));
spikeTrain(spike_index)=1;
figure;
plot(linspace(0,maxtime,length(data)),spikeTrain)
title('spike train')
xlabel('time(s)')
ylabel('spike')
%% Estimated firing rate
% Rectangular window
deltaT=0.5;
windowlength=floor(deltaT*par.sr);
window=ones(1,windowlength);
stimulatedFR=conv(spikeTrain,window);
% plot the resualt
figure;
plot(linspace(0,length(stimulatedFR)/par.sr,length(stimulatedFR)),stimulatedFR);
title('Estimated firing rate with rectangular window')
xlabel('time(s)')
ylabel('firing rate')
%% Gaussian window
delta_t=0.7;
windowlengthG=delta_t*par.sr;
alpha=3;%Standard Daviation(SD)
windowG=gausswin(windowlengthG,alpha);
stimulatedFRG=conv(spikeTrain,windowG);
figure;
plot(linspace(0,length(stimulatedFRG)/par.sr,length(stimulatedFRG)),stimulatedFRG);
title('Estimated firing rate with Gaussian window')
xlabel('time(s)')
ylabel('firing rate')
%%
b=spike_index/par.sr;
b = sort(b);
%% Inter Spike Interval distribution & Coefficient of Variation (CV)
ISI=diff(b);
histogram(ISI)
title('Inter Spike Interval distribution')
ylabel('Inter Spike Interval distribution')
mean_ISI=mean(b);
std_ISI=std(b);
CV=(std_ISI/mean_ISI)*100;
fprintf('Coefficient of Variation (CV):%.2f%%\n',CV');
%% the spike-train auto correlation function

% Compute auto-correlation using xcorr
lag = -(length(b)-1):(length(b)-1);
auto_corr = xcorr(b,b, 'coeff');

% Plot the auto-correlation function
stem(lag, auto_corr);
xlabel('Lag');
ylabel('Auto-correlation');
title('Auto-correlation of Spike Train');
%% Fano Factor
bin_size=0.5;
bin_edge=0:bin_size:max(b);
spike_counts=histc(b,bin_edge);
mean_spike_count=mean(spike_counts);
var_spike_count=var(spike_counts);
fano_factor=var_spike_count/mean_spike_count;
disp(['Fano Factor:',num2str(fano_factor)]);
%%
clear all
close all
clc
%% load data
load 'times_090425blk10_ch66.mat'
%% multi trial analysis
par.scales=4;
data=spikes(:,:)*par.scales;
figure;
plot(data)
title('Raw data');
%% fit sampling rate in PSD
pwelch(data,[],[],[],par.sr);
%% filtering
order=3;
cutoff_freq=1100;%Hz
[b,a] = butter(order,cutoff_freq/(par.sr/2),'high');
filtered_data=filtfilt(b,a,data);
figure;
plot(filtered_data)
title('filtered data')
hold on
% spike detection
sig=zeros(1,length(filtered_data(1,:)));
for h=1:length(filtered_data(1,:))
    sig(1,h)=median(abs(filtered_data(1,h))/0.6745);
end
plot(ones(length(filtered_data(:,:)),1)*5*sig,'LineWidth',2)
%%
spike_train=zeros(length(filtered_data(:,1)),length(filtered_data(1,:)));
for f=1:length(filtered_data(1,:))
    spike_train(:,f)=(filtered_data(:,f))>5*sig(f);
end
%%
spike=zeros(length(spike_train(:,1)),length(spike_train(1,:)));
spike_index=zeros(length(spike_train(:,1)),length(spike_train(1,:)));
for k=1:length(spike_train(1,:))
    oneind = find(spike_train(:,k));
    for n = 1:size(oneind,1)-1
        if (oneind(n+1) == oneind(n) + 1)
            spike_train(oneind(n+1),k) = 0;
        end
    end
    spike_train_New(:,k) = spike_train(:,k) + (k - 1);
    plot(spike_train_New(:,k))
    hold on
end

for v = 1:size(spike_train,1)
    if  mean(spike_train(v,:)) >= 0.4
        spike_train_mean(v,1) = 1;
    else
        spike_train_mean(v,1) = 0;
    end
end
%% created spike train
maxtime=length(data)/par.sr;
[spike,spike_index]=findpeaks(spike_train_mean,'MinPeakDistance',300);
spikeTrain=zeros(1,length(spike_train_mean));
spikeTrain(spike_index)=1;
figure;
plot(linspace(0,maxtime,length(data)),spikeTrain)
title('spike train')
xlabel('time(s)')
ylabel('spike')
%% created spike train
spikeTrain=zeros(1,length(data));
spikeTrain(spike_index)=1;
figure;
stem(linspace(0,maxtime,length(data)),spikeTrain)
title('spike train')
xlabel('time(s)')
ylabel('spike')
%% Estimated firing rate
% Rectangular window
deltaT=1;
windowlength=floor(deltaT*par.sr);
window=ones(1,windowlength);
stimulatedFR=conv(spikeTrain,window);
% plot the resualt
figure;
plot(linspace(0,length(stimulatedFR)/par.sr,length(stimulatedFR)),stimulatedFR);
title('Estimated firing rate with rectangular window')
xlabel('time(s)')
ylabel('firing rate')
%% Gaussian window
delta_t=1;
windowlengthG=delta_t*par.sr;
alpha=6;%Standard Daviation(SD)
windowG=gausswin(windowlengthG,alpha);
stimulatedFRG=conv(spikeTrain,windowG);
figure;
plot(linspace(0,length(stimulatedFRG)/par.sr,length(stimulatedFRG)),stimulatedFRG);
title('Estimated firing rate with Gaussian window')
xlabel('time(s)')
ylabel('firing rate')
%%
b=spike_index/par.sr;
b = sort(b);
%% Inter Spike Interval distribution & Coefficient of Variation (CV)
ISI=diff(b);
histogram(ISI)
title('Inter Spike Interval distribution')
ylabel('Inter Spike Interval distribution')
mean_ISI=mean(b);
std_ISI=std(b);
CV=(std_ISI/mean_ISI)*100;
fprintf('Coefficient of Variation (CV):%.2f%%\n',CV');
%% the spike-train auto correlation function

% Compute auto-correlation using xcorr
lag = -(length(b)-1):(length(b)-1);
auto_corr = xcorr(b,b, 'coeff');

% Plot the auto-correlation function
stem(lag, auto_corr);
xlabel('Lag');
ylabel('Auto-correlation');
title('Auto-correlation of Spike Train');
%% Fano Factor
bin_size=0.5;
bin_edge=0:bin_size:max(b);
spike_counts=histc(b,bin_edge);
mean_spike_count=mean(spike_counts);
var_spike_count=var(spike_counts);
fano_factor=var_spike_count/mean_spike_count;
disp(['Fano Factor:',num2str(fano_factor)]);
