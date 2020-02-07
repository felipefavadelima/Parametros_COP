%Código que gera sinal com SLP_PSD conhecido
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

%Eixo do Tempo
dt = 1/Fs;
t = 0:dt:T-dt;

%Gera sinal com
SLP_PSD_LF = -2;
SLP_PSD_HF = -5;

%LF - PSD em Baixas frequências (0.05Hz até 1Hz)
logxl = log10(0.05):0.001:log10(1);
logyl = SLP_PSD_LF.*logxl;
xl = 10.^logxl;
yl = 10.^logyl;

%HF -PSD em Altas frequências (1Hz até 5)
logxh = log10(1+0.001):0.001:log10(5);
logyh = SLP_PSD_HF.*logxh;
xh = 10.^logxh;
yh = 10.^logyh;

%PSD completo
flog = [0 xl xh Fs/2];
PSD = [0.0001 yl yh 0.0001];
PSD = interp1(flog,PSD,f);

%Plot logxlog do PSD simulado
subplot(2,1,1);
plot(f,log10(PSD));
xlabel('Freq(Hz)');
ylabel('log(PSD)');
vline(0.05,'b','0.05Hz');
vline(1,'b','1Hz');
vline(5,'b','5Hz');
xlim([0 10]);
subplot(2,1,2);
plot(log10(f),log10(PSD));
xlabel('log(Freq(Hz))');
ylabel('log(PSD)');
vline(log10(0.05),'b','0.05Hz');
vline(log10(1),'b','1Hz');
vline(log10(5),'b','5Hz');

%Gera sinal com PSD especificado
FFTabs = sqrt(PSD);
FFTang = pi*(randn(1,length(FFTabs)));
FFTabs = [FFTabs fliplr(FFTabs(2:end))];
FFTabs = FFTabs';
FFTang = [FFTang -fliplr(FFTang(2:end))];
FFTang = FFTang';
[FFTr,FFTi] = pol2cart(FFTang,FFTabs);
Sinal_simulado = ifft(FFTr+1i*FFTi);
Sinal_simulado = real(Sinal_simulado);

%Sinal simulado
Sinal_simulado(end+1) = Sinal_simulado(end);

%Plot fftSinal_simulado log x log
logfftSinal_simulado = log10(abs(fft(Sinal_simulado)).^2);
logfftSinal_simulado = logfftSinal_simulado(1:end/2);
plot(log10(f),logfftSinal_simulado);
vline(log10(0.05),'b','0.05Hz');
vline(log10(1),'b','1Hz');
vline(log10(5),'b','5Hz');
ylabel('FFT\{Sinal_simulado \}^2');
xlabel('frequência(Hz)');