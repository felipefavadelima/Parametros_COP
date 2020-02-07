%Código que gera sinal com f50p_PSD conhecido
%por: Felipe Fava de Lima @2020
%email: felipefavadelima@gmail.com
clear all
close all

%Parâmetros do sinal simulado
T = 90;
Fs = 100;
N = Fs*T;

%Eixo das Frequências
df = 1/T;
f = 0:df:(Fs/2)-df;

%Eixos do Tempo
dt = 1/Fs;
t = 0:dt:T-dt;

%Gera sinal com f50p = 2.5725Hz
%PSD com formato triangular de 0.15Hz a 5Hz
f50p_PSD_APver = ((5+0.15)/2);
Idxs = intersect(find(f>0.15),find(f<5));
FFTabs = 0:(1/(length(Idxs)/2)):1;
FFTabs = [FFTabs fliplr(FFTabs)];
FFTabs = [zeros(Idxs(1)-2,1); FFTabs'; zeros(length(f)-Idxs(end),1)];
FFTang = pi*(randn(length(FFTabs),1));
FFTabs = [FFTabs' fliplr(FFTabs(2:end)')];
FFTabs = sqrt(FFTabs');
FFTang = [FFTang' -fliplr(FFTang(2:end)')];
FFTang = FFTang';
[FFTr,FFTi] = pol2cart(FFTang,FFTabs);
Sinal_simulado = ifft([FFTr+i*FFTi]);

%Sinal simulado
Sinal_simulado  = 1000*real(Sinal_simulado (1:N));

%Plot do sinal simulado
figure
subplot(2,1,1);
plot(t,Sinal_simulado );
ylabel('Sinal_simulado ');
xlabel('Tempo(s)');
subplot(2,1,2);
fftAP = abs(fft(Sinal_simulado ).^2);
plot(f,fftAP(1:length(f)));
ylabel('FFT\{Sinal_simulado \}^2');
xlabel('frequência(Hz)');