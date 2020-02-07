%Código para o cálculo dos parâmetros f50p_PSD e SLP_PSD
%por: Felipe Fava de Lima @2020
%email: felipefavadelima@gmail.com
clear all
close all
%Entradas
%   AP - Vetor de dados da localização do COP na direção AP
%   ML - Vetor de dados da localização do COP na direção ML

%Parâmetros de saída
%   f50p_PSD_AP
%   f50p_PSD_AP
%   SLP_PSD_AP_HF
%   SLP_PSD_AP_LF
%   SLP_PSD_ML_HF
%   SLP_PSD_ML_LF

%%
%Resample
%Filtro passa-baixas Butterworth 10Hz fase nula 4ªordem
%[b,a] = butter(4,(10/(1000/2)));
%AP = filtfilt(b,a,AP);
%ML = filtfilt(b,a,ML);

%%
%Re-amostragem do sinal de fs= 1000Hz para fs= 100Hz
%AP = resample(AP,100,1000);
%ML = resample(ML,100,1000);

%%
%Frequência de amostragem em Hz
Fs = 100;

%Número de amostras adquiridas
N = length(AP);

%Duração da aquisição em segundos
T = N * (1/Fs);

%%
%Equações (1) e (2) de (PRIETO et al.,  1996)
ML0 = ML - mean(ML);
AP0 = AP - mean(AP);

%Estimador de PSD Pwelch
WinT_s = 20;
WinSize = WinT_s*Fs;
Psd_df =  1/WinT_s;
Win = hann(WinSize);
[Psd_AP Psd_freqs] = pwelch(AP0,Win,WinSize/2,WinSize,Fs);
Psd_ML = pwelch(ML0,Win,WinSize/2,WinSize,Fs);

%%
%PARÂMETROS f50p_PSD (0.15 Hz a 5 Hz)
%AP
cum_Area_Psd_AP = cumsum(PSD_AP.*Psd_df);
freq_Idxs = intersect(find(PSD_freqs < 5), find(PSD_freqs > 0.15));
Area_Psd_AP = cum_Area_Psd_AP(freq_Idxs(end)) - ...
    cum_Area_Psd_AP(freq_Idxs(1));
half_Area_Psd_AP = 0.5 * Area_Psd_AP;
Idx_half_Area_Psd_AP = find(cum_Area_Psd_AP < half_Area_Psd_AP +...
    cum_Area_Psd_AP(freq_Idxs(1)));
Idx_half_Area_Psd_AP = Idx_half_Area_Psd_AP(end)+1;
f50p_PSD_AP = PSD_freqs(Idx_half_Area_Psd_AP(end));
%ML
cum_Area_Psd_ML = cumsum(PSD_ML.*Psd_df);
freq_Idxs = intersect(find(PSD_freqs < 5), find(PSD_freqs > 0.15));
Area_Psd_ML = cum_Area_Psd_ML(freq_Idxs(end)) - ...
    cum_Area_Psd_ML(freq_Idxs(1));
half_Area_Psd_ML = 0.5 * Area_Psd_ML;
Idx_half_Area_Psd_ML = find(cum_Area_Psd_ML < half_Area_Psd_ML +...
    cum_Area_Psd_ML(freq_Idxs(1)));
Idx_half_Area_Psd_ML = Idx_half_Area_Psd_ML(end)+1;
f50p_PSD_ML = PSD_freqs(Idx_half_Area_Psd_ML(end));

%%
%PARÂMETROS SLP_PSD

%Baixas frequências (0.05Hz até 1Hz)
Psd_low_Freqs_Idxs =  intersect(find(Psd_freqs<=1),...
    find(Psd_freqs>=0.05));
Psd_low_Freqs = Psd_freqs(Psd_low_Freqs_Idxs);
Psd_AP_low_freqs = Psd_AP(Psd_low_Freqs_Idxs);
Psd_ML_low_freqs = Psd_ML(Psd_low_Freqs_Idxs);

%SLP_PSD_AP_LF
SLP_PSD_AP_LF = polyfit(log10(Psd_low_Freqs),...
    log10(Psd_AP_low_freqs),1);
SLP_PSD_AP_LF = SLP_PSD_AP_LF(1);

%SLP_PSD_ML_LF
SLP_PSD_ML_LF = polyfit(log10(Psd_low_Freqs),...
    log10(Psd_ML_low_freqs),1);
SLP_PSD_ML_LF = SLP_PSD_ML_LF(1);

%Altas frequências (1Hz até 5)
Psd_high_Freqs_Idxs =  intersect(find(Psd_freqs<=5),...
    find(Psd_freqs>=1));
Psd_high_Freqs = Psd_freqs(Psd_high_Freqs_Idxs);
Psd_AP_high_freqs = Psd_AP(Psd_high_Freqs_Idxs);
Psd_ML_high_freqs = Psd_ML(Psd_high_Freqs_Idxs);

%SLP_PSD_AP_HF
SLP_PSD_AP_HF = polyfit(log10(Psd_high_Freqs),...
    log10(Psd_AP_high_freqs),1);
SLP_PSD_AP_HF = SLP_PSD_AP_HF(1);

%SLP_PSD_ML_HF
SLP_PSD_ML_HF = polyfit(log10(Psd_high_Freqs),...
    log10(Psd_ML_high_freqs),1);
SLP_PSD_ML_HF = SLP_PSD_ML_HF(1);