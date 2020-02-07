%Código que gera COP com AREACE E PHI conhecidos
%por: Felipe Fava de Lima @2020
%email: felipefavadelima@gmail.com
clear all
close all

%%
%Parâmetros do sinal
T = 90;
fs = 1000;
N = fs*T;

%Parâmetros da elipse 95%
AP_std = 2;  
ML_std = 10;

%Phi simulado
PHI_simulado = 90;

%AREACE da elipse 95% simulado
AREACE_simulado = pi*(3*AP_std)*(3*ML_std);

%%
%Gera AP ML por variável aleatória com
%distribuição normal com std() AP_std e ML_std
AP_0 = AP_std* randn(N,1);
ML_0 = ML_std*randn(N,1);

%Rotaciona o sinal do COP
if(AP_std > ML_std)
    alpha = ((90-PHI_simulado)*pi)/180;
else
    alpha = ((180-PHI_simulado)*pi)/180;
end
%Sinal do COP simulado
AP = AP_0.*cos(alpha) - ML_0.*sin(alpha);
ML = AP_0.*sin(alpha) + ML_0.*cos(alpha);

%%
%Plot do COP simulado
plot(ML,AP);
dAP = abs(max(AP) - min(AP));
dML = abs(max(ML) - min(ML));
dmax = max([dAP dML]);
daxis = 1.2 * (dmax/2);
xlim([mean(ML)-daxis mean(ML)+daxis]);
ylim([mean(AP)-daxis mean(AP)+daxis]);
hold on

%Plot da elipse 95%
w = 0:(2*pi)/1000:2*pi;
EliML0 = 3*ML_std*cos(w);
EliAP0 = 3*AP_std*sin(w);
EliAP = EliAP0.*cos(alpha) - EliML0.*sin(alpha);
EliML = EliAP0.*sin(alpha) + EliML0.*cos(alpha);
plot(EliML,EliAP);
hline(0,'r','ML');
vline(0,'k','AP');