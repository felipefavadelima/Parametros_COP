%Código para o cálculo dos parâmetros CI da Multi scale entropy
%por: Felipe Fava de Lima @2020
%email: felipefavadelima@gmail.com
clear all
close all
%Entradas
%   AP - Vetor de dados da localização do COP na direção AP
%   ML - Vetor de dados da localização do COP na direção ML

%Parâmetros de saída
%   CI_AP
%   CI_ML

%%
%Pré-processamento
%Filtro passa-faixas Butterworth 0.05Hz a 20Hz fase nula 4ªordem
%[b,a] = butter(4,([0.05 20]/(1000/2)));
%AP = filtfilt(b,a,AP);
%ML = filtfilt(b,a,ML);

%%
%Resample
%Re-amostragem do sinal de fs= 1000Hz para fs= 100Hz
%AP = resample(AP,100,1000);
%ML = resample(ML,100,1000);

%%
%Frequência de amostragem em Hz
Fs = 100;

%Número de amostras adquiridas
S = length(AP);

%Duração da aquisição em segundos
T = S * (1/Fs);

%%
%Parâmetros
m=2;
r = 0.2;
MinScale_s = 0.03;
MaxScale_s = 0.33;


%%
%Escalas
NMinScale = MinScale_s*Fs;
NMaxScale = MaxScale_s*Fs;
dScales = round((NMaxScale-NMinScale)/10);
Scales = NMinScale:dScales:NMaxScale;

%%
%AP
Signal = AP;
SampEn = zeros(length(Scales),1);
for ScaleIdx=1:length(Scales)
    Scale = Scales(ScaleIdx);
    %Gera sinal y na escala de tempo
    Signal = Signal(1:Scale*fix(length(Signal)/Scale));
    u = reshape(Signal,Scale,fix(length(Signal)/Scale));
    u = mean(u,1);
    %Calcula Threshold
    Thrs=r*std(u);
    S = length(u);
    %%
    %Construção vetores x(i,m) = [y(i),...,y(i+m+1)];
    x = [];
    for delay=1:m+1
       Temp = u(delay:1:S-(m-delay)-1);
       x=[x;Temp];
    end
    %%
    %Cálculo de 
    %   d[x(i),x(j)] = max(x(i) - x(j));
    %   phis  = estimativa de ocorrências de (d[x(i),x(j)]<Thrs)
    %   PHI = média de phis
    %para m
    xj = x(1:end-1,:);
    for i=1:length(xj)
        xi = repmat(xj(:,i),1,length(xj));
        d = max(abs(xi-xj));
        Nmatchs = sum(d<Thrs);
        phis_m(i) = (Nmatchs-1)/(length(xj)-1);
    end
    Phi_m = mean(phis_m);
    %para m+1
    xj = x;
    for i=1:length(xj)
        xi = repmat(xj(:,i),1,length(xj));
        d = max(abs(xi-xj));
        Nmatchs = sum(d<Thrs);
        phis_mp1(i) = (Nmatchs-1)/(length(xj)-1);
    end
    Phi_mp1 = mean(phis_mp1);
    SampEn(ScaleIdx) = -log(Phi_1/Phi_m);
end
CI_AP = sum(SampEn);

%%
%ML
Signal = ML;
SampEn = zeros(length(Scales),1);
for ScaleIdx=1:length(Scales)
    Scale = Scales(ScaleIdx);
    %Gera sinal y na escala de tempo
    Signal = Signal(1:Scale*fix(length(Signal)/Scale));
    u = reshape(Signal,Scale,fix(length(Signal)/Scale));
    u = mean(u,1);
    %Calcula Threshold
    Thrs=r*std(u);
    S = length(u);
    %%
    %Construção vetores x(i,m) = [y(i),...,y(i+m+1)];
    x = [];
    for delay=1:m+1
       Temp = u(delay:1:S-(m-delay)-1);
       x=[x;Temp];
    end
    %%
    %Cálculo de 
    %   d[x(i),x(j)] = max(x(i) - x(j));
    %   phis  = estimativa de ocorrências de (d[x(i),x(j)]<Thrs)
    %   PHI = média de phis
    %para m
    xj = x(1:end-1,:);
    for i=1:length(xj)
        xi = repmat(xj(:,i),1,length(xj));
        d = max(abs(xi-xj));
        Nmatchs = sum(d<Thrs);
        phis_m(i) = (Nmatchs-1)/(length(xj)-1);
    end
    Phi_m = mean(phis_m);
    %para m+1
    xj = x;
    for i=1:length(xj)
        xi = repmat(xj(:,i),1,length(xj));
        d = max(abs(xi-xj));
        Nmatchs = sum(d<Thrs);
        phis_mp1(i) = (Nmatchs-1)/(length(xj)-1);
    end
    Phi_mp1 = mean(phis_mp1);
    SampEn(ScaleIdx) = -log(Phi_mp1/Phi_m);
end
CI_ML = sum(SampEn);