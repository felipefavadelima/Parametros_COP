%Código para o cálculo dos parâmetros RDIST, MVELO, AREACE e PHI
%por: Felipe Fava de Lima @2020
%email: felipefavadelima@gmail.com
clear all
close all

%Entradas
%   AP[n] - Vetor de dados da localização do COP na direção AP
%   ML[n] - Vetor de dados da localização do COP na direção ML

%Parâmetros de saída
%   RDIST_AP 
%   RDIST_ML
%   RDIST_AP_LF 
%   RDIST_ML_LF 
%   RDIST_AP_HF 
%   RDIST_ML_HF
%   MVELO_AP
%   MVELO_ML 
%   MVELO_AP_LF 
%   MVELO_ML_LF 
%   MVELO_AP_HF 
%   MVELO_ML_HF 
%   AREACE 
%   AREACE_HF
%   AREACE_LF
%   PHI
%   PHI_LF
%   PHI_HF

%%
%Re-amostragem do sinal de fs= 1000Hz para fs= 100Hz
AP = resample(AP,100,1000);
ML = resample(ML,100,1000);

%%
%Frequência de amostragem em Hz
Fs = 100;

%Número de amostras adquiridas
N = length(AP);

%Duração da aquisição em segundos
T = N * (1/Fs);

%%
%Filtro passa-baixas Butterworth 10Hz fase nula 4ªordem
[b,a] = butter(4,(10/(Fs/2)));
AP = filtfilt(b,a,AP);
ML = filtfilt(b,a,ML);

%Filtro passa-baixas Butterworth 0.5Hz fase nula 4ªordem
[b,a] = butter(4,(0.5/(Fs/2)));
AP_LF = filtfilt(b,a,AP);
ML_LF = filtfilt(b,a,ML);

%Filtro passa-altas Butterworth 0.5Hz fase nula 4ªordem
[b,a] = butter(4,(0.5/(Fs/2)),'high');
AP_HF = filtfilt(b,a,AP);
ML_HF = filtfilt(b,a,ML);

%%
%Equações (1) e (2) de (PRIETO et al.,  1996)
ML0 = ML - mean(ML);
AP0 = AP - mean(AP);
ML0_LF = ML_LF - mean(ML_LF);
AP0_LF = AP_LF - mean(AP_LF);
ML0_HF = ML_HF - mean(ML_HF);
AP0_HF = AP_HF - mean(AP_HF);

%%
%Parâmetros RDIST
%Equação (6) de (PRIETO et al.,  1996)
RDIST_AP = rms(AP0);
RDIST_ML = rms(ML0);
RDIST_AP_LF = rms(AP0_LF);
RDIST_ML_LF = rms(ML0_LF);
RDIST_AP_HF = rms(AP0_HF);
RDIST_ML_HF = rms(ML0_HF);

%%
%Distância total percorrida pelo COP
%Equações (8) e (9) de (PRIETO et al.,  1996)
TOTEX_AP = sum(diff(AP0));
TOTEX_ML = sum(diff(ML0));
TOTEX_AP_LF = sum(diff(AP0_LF));
TOTEX_ML_LF = sum(diff(ML0_LF));
TOTEX_AP_HF = sum(diff(AP0_HF));
TOTEX_ML_HF = sum(diff(ML0_HF));

%%
%Parâmetros MVELO
%Equações (10) e (11) de (PRIETO et al.,  1996)
MVELOAP = TOTEX_AP/T;
MVELOML = TOTEX_ML/T;
MVELOAP_LF = TOTEX_AP_LF/T;
MVELOML_LF = TOTEX_ML_LF/T;
MVELOAP_HF = TOTEX_AP_HF/T;
MVELOML_HF = TOTEX_ML_HF/T;

%%
%Parâmetros AREACE e PHI
%Variâncias e covariância
VAR_AP = sum(AP0.^2)/(N-1);
VAR_ML = sum(ML0.^2)/(N-1);
COV_APML = sum(AP0.*ML0)/(N-1);

VAR_AP_LF = sum(AP0_LF.^2)/(N-1);
VAR_ML_LF = sum(ML0_LF.^2)/(N-1);
COV_APML_LF = sum(AP0_LF.*ML0_LF)/(N-1);

VAR_AP_HF = sum(AP0_HF.^2)/(N-1);
VAR_ML_HF = sum(ML0_HF.^2)/(N-1);
COV_APML_HF = sum(AP0_HF.*ML0_HF)/(N-1);

%AREA CE
%Equações (16), (17) e (18) de (PRIETO et al.,  1996)
%Obs:Aparentemente a equação (18) de (PRIETO et al.,  1996)
%deve ser corrigida para: 
%AREA_CE = (F05[2,n-2])^2*pi*sqrt(SAP^2SML^2 - SAPML^2)
AREACE = 9*pi*sqrt(VAR_AP*VAR_ML - COV_APML^2);
AREACE_LF = 9*pi*sqrt(VAR_AP_LF*VAR_ML_LF - COV_APML_LF^2);
AREACE_HF = 9*pi*sqrt(VAR_AP_HF*VAR_ML_HF - COV_APML_HF^2);

%%
%PHI
%Procedimento baseado em (KUO et al., 1998)
%Matriz covariância
C = [[VAR_ML COV_APML];[COV_APML VAR_AP]]; 
C_LF = [[VAR_ML_LF COV_APML_LF];[COV_APML_LF VAR_AP_LF]]; 
C_HF = [[VAR_ML_HF COV_APML_HF];[COV_APML_HF VAR_AP_HF]]; 

MLv = [1;0];                                %Eixo ML
                                
                                            %PHI
[V,e] = eig(C)                              %Calcula autovetor 
                                            %e autovalor
if(e(1,1) > e(2,2))                         %Escolhe maior autovetor
    ElipAx = V(:,1);
else
    ElipAx = V(:,2);
end
if(ElipAx(2) < 0 && ElipAx(1) < 0)          %Se no terceiro quadrante
    ElipAx = - ElipAx;                      %Envia para primeiro
                                            %quadrante
elseif(ElipAx(1) > 0 && ElipAx(2) < 0)      %Se no quarto quadrante
    ElipAx = - ElipAx;                      %Envia para primeiro
end
cosPHI = sum(ElipAx.*MLv)/...
    (norm(ElipAx) * norm(MLv));             %calcula cos do ângulo entre
PHI = acos(cosPHI);                         %maior eixo da elipse e ML
PHI = (PHI*180)/pi;                         %Tranforma para graus

                                            %PHI_LF
[V_LF,e_LF] = eig(C_LF)                     %Calcula autovetor e
                                            %autovalor
if(e_LF(1,1) > e_LF(2,2))                   %Escolhe maior autovetor
    ElipAx_LF = V_LF(:,1);
else
    ElipAx_LF = V_LF(:,2);
end
if(ElipAx_LF(2) < 0 && ElipAx_LF(1) < 0)    %Se no terceiro quadrante
    ElipAx_LF = - ElipAx_LF;                %Envia para primeiro
                                            %quadrante
elseif(ElipAx_LF(1) > 0 && ElipAx_LF(2) < 0)%Se no quarto quadrante
    ElipAx_LF = - ElipAx_LF;                %Envia para primeiro
end
cosPHI_LF = sum(ElipAx_LF.*MLv)/...
    (norm(ElipAx_LF) * norm(MLv));          %calcula cos do ângulo entre
PHI_LF = acos(cosPHI_LF);                   %maior eixo da elipse e ML
PHI_LF = (PHI_LF*180)/pi;                   %Tranforma para graus

                                            %PHI_HF
[V_HF,e_HF] = eig(C_HF)                     %Calcula autovetor e
                                            %autovalor
if(e_HF(1,1) > e_HF(2,2))                   %Escolhe maior autovetor
    ElipAx_HF = V_HF(:,1);
else
    ElipAx_HF = V_HF(:,2);
end
if(ElipAx_HF(2) < 0 && ElipAx_HF(1) < 0)    %Se no terceiro quadrante
    ElipAx_HF = - ElipAx_HF;                %Envia para primeiro
                                            %quadrante
elseif(ElipAx_HF(1) > 0 && ElipAx_HF(2) < 0)%Se no quarto quadrante
    ElipAx_HF = - ElipAx_HF;                %Envia para primeiro
end
cosPHI_HF = sum(ElipAx_HF.*MLv)/...
    (norm(ElipAx_HF) * norm(MLv));          %calcula cos do ângulo entre
PHI_HF = acos(cosPHI_HF);                   %maior eixo da elipse e ML
PHI_HF = (PHI_HF*180)/pi;                   %Tranforma para graus