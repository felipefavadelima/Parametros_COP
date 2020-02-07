%Código que gera sinal com características de entropia conhecidas
%por: Felipe Fava de Lima @2020
%email: felipefavadelima@gmail.com
clear all
close all

%Parâmetros do sinal
T = 90;
Fs = 100;
N = Fs*T
dt = 1/Fs;
t = 0:dt:T-dt;

%Processo MIX (Joshua,2000)   
Sinal_simulado = 3*sin(2*pi*fsin*t);
P = 0;
IdxMix = randperm(N,P*N);
RandN = randn(P*N,1);
Sinal_simulado(IdxMix) = RandN;

%Plot Sinal_simulado
plot(t,Sinal_simulado);