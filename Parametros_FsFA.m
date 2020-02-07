%Código para o cálculo dos parâmetros alpha_fsFA
%por: Felipe Fava de Lima @2020
%email: felipefavadelima@gmail.com
clear all
close all

%Entradas
%   AP - Vetor de dados da localização do COP na direção AP
%   ML - Vetor de dados da localização do COP na direção ML

%Parâmetros de saída
%   ?_FSFA_AP_LF
%   ?_FSFA_AP_HF
%   ?_FSFA_ML_LF
%   ?_FSFA_ML_HF
%   ?_FSFA_vAP_LF
%   ?_FSFA_vAP_HF
%   ?_FSFA_vML_LF
%   ?_FSFA_vML_HF

%%
%Filtro passa-baixas Butterworth 20Hz fase nula 4ªordem
%[b,a] = butter(4,(20/(1000/2)));
%AP = filtfilt(b,a,AP);
%ML = filtfilt(b,a,ML);

%%
%Re-amostragem do sinal de fs= 1000Hz para fs= 40Hz
%AP = resample(AP,40,1000);
%ML = resample(ML,40,1000);

%Frequência de amostragem em Hz
Fs = 40;

%Velocidade de AP e ML
vAP = diff(AP)*Fs;
vML = diff(ML)*Fs;

%%
%Integral dos sinais
APint = cumsum(AP-mean(AP));
MLint = cumsum(ML-mean(ML));
vAPint = cumsum(vAP-mean(vAP));
vMLint = cumsum(vML-mean(vML));

%%
%Eixo pseudo exponencial de tamanhos janelas
npts = 10;
N = min([length(APint) length(vAPint)]);
Wsizesl = round(logspace(log10(4),log10(20),npts));
Wsizesh = round(logspace(log10(80),log10(N/4),npts));
Wsizes = [Wsizesl Wsizesh];                             %Tamando das 
                                                        %janelas
WIdxsWl = 1:npts;                                       %Índices LF
WIdxsWh = npts+1:2*npts;                                %Índices HF

%%
%AP
Signal = APint;
N = length(Signal);
F = zeros(length(Wsizes),1);
for nWsize=1:length(Wsizes)
    W = Wsizes(nWsize);
    Wins = reshape(Signal(1:W*fix(N/W)),[W,fix(N/W)]);  %Janelamento
    RmsDetWins = zeros(fix(N/W),1);
    for nWin=1:fix(N/W)
        Win = detrend(Wins(:,nWin));                    %Detrend das 
                                                        %janelas
        RmsDetWins(nWin) = rms(Win);                    %Rms das janelas
    end
    F(nWsize) = rms(RmsDetWins);                        %RMS do sinal
end

%?_FSFA_AP_LF
ALPHA_FSFA_AP_LF = polyfit(log10(Wsizes(WIdxsWl)'),...  %Polyfit
    log10(F(WIdxsWl)),1);
ALPHA_FSFA_AP_LF = ALPHA_FSFA_AP_LF(1);

%?_FSFA_AP_HF
ALPHA_FSFA_AP_HF = polyfit(log10(Wsizes(WIdxsWh)'),...  %Polyfit
    log10(F(WIdxsWh)),1);
ALPHA_FSFA_AP_HF = ALPHA_FSFA_AP_HF(1);

%%
%ML
Signal = MLint;
N = length(Signal);
F = zeros(length(Wsizes),1);
for nWsize=1:length(Wsizes)
    W = Wsizes(nWsize);
    Wins = reshape(Signal(1:W*fix(N/W)),[W,fix(N/W)]);  %Janelamento
    RmsDetWins = zeros(fix(N/W),1);
    for nWin=1:fix(N/W)
        Win = detrend(Wins(:,nWin));                    %Detrend das 
                                                        %janelas
        RmsDetWins(nWin) = rms(Win);                    %Rms das janelas
    end
    F(nWsize) = rms(RmsDetWins);                        %RMS do sinal
end

%?_FSFA_ML_LF
ALPHA_FSFA_ML_LF = polyfit(log10(Wsizes(WIdxsWl)'),...  %Polyfit
    log10(F(WIdxsWl)),1);
ALPHA_FSFA_ML_LF = ALPHA_FSFA_ML_LF(1);

%?_FSFA_ML_HF
ALPHA_FSFA_ML_HF = polyfit(log10(Wsizes(WIdxsWh)'),...  %Polyfit
    log10(F(WIdxsWh)),1);
ALPHA_FSFA_ML_HF = ALPHA_FSFA_ML_HF(1);

%%
%vAP
Signal = vAPint;
N = length(Signal);
F = zeros(length(Wsizes),1);
for nWsize=1:length(Wsizes)
    W = Wsizes(nWsize);
    Wins = reshape(Signal(1:W*fix(N/W)),[W,fix(N/W)]);  %Janelamento
    RmsDetWins = zeros(fix(N/W),1);
    for nWin=1:fix(N/W)
        Win = detrend(Wins(:,nWin));                    %Detrend das 
                                                        %janelas
        RmsDetWins(nWin) = rms(Win);                    %Rms das janelas
    end
    F(nWsize) = rms(RmsDetWins);                        %RMS do sinal
end

%?_FSFA_vAP_LF
ALPHA_FSFA_vAP_LF = polyfit(log10(Wsizes(WIdxsWl)'),...  %Polyfit
    log10(F(WIdxsWl)),1);
ALPHA_FSFA_vAP_LF = ALPHA_FSFA_vAP_LF(1);

%?_FSFA_vAP_HF
ALPHA_FSFA_vAP_HF = polyfit(log10(Wsizes(WIdxsWh)'),...  %Polyfit
    log10(F(WIdxsWh)),1);
ALPHA_FSFA_vAP_HF = ALPHA_FSFA_vAP_HF(1);

%%
%vML
Signal = vMLint;
N = length(Signal);
F = zeros(length(Wsizes),1);
for nWsize=1:length(Wsizes)
    W = Wsizes(nWsize);
    Wins = reshape(Signal(1:W*fix(N/W)),[W,fix(N/W)]);  %Janelamento
    RmsDetWins = zeros(fix(N/W),1);
    for nWin=1:fix(N/W)
        Win = detrend(Wins(:,nWin));                    %Detrend das 
                                                        %janelas
        RmsDetWins(nWin) = rms(Win);                    %Rms das janelas
    end
    F(nWsize) = rms(RmsDetWins);                        %RMS do sinal
end

%?_FSFA_vML_LF
ALPHA_FSFA_vML_LF = polyfit(log10(Wsizes(WIdxsWl)'),...  %Polyfit
    log10(F(WIdxsWl)),1);
ALPHA_FSFA_vML_LF = ALPHA_FSFA_vML_LF(1);

%?_FSFA_vML_HF
ALPHA_FSFA_vML_HF = polyfit(log10(Wsizes(WIdxsWh)'),...  %Polyfit
    log10(F(WIdxsWh)),1);
ALPHA_FSFA_vML_HF = ALPHA_FSFA_vML_HF(1);