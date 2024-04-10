clear all; close all; clc
%A1
T=10^-3;
over=10;
Ts=T/over;
A=4;
a = [0,0.5,1];

phit_matrix = [];

figure;
 hold on;
 grid on;
for i=1:length(a)
   [phi,t] = srrc_pulse(T,over,A,a(i));
   phit_matrix = [phit_matrix; phi];
   plot(t,phi);
end
hold off;
title('SRRC Pulses');
xlabel('Time');
ylabel('Amplitude');
legend('a=0','a=0.5','a=1');

%A2
%a)
N=2048;
Fs = 1/Ts;
f_axis = (-0.5:1/N:0.5-1/N);
F_axis = Fs*f_axis;
PHIF_matrix = [];
PHI_psd_matrix=[];
figure;
hold on;
grid on;
for i=1:height(phit_matrix)
    PHI = fftshift(fft(phit_matrix(i,:),N))*Ts;
    PHIF_matrix =[PHIF_matrix; PHI];
    PHI_psd = abs(PHI).^2;
    PHI_psd_matrix=[PHI_psd_matrix;PHI_psd];
    plot(F_axis,PHI_psd);
end
hold off;
title("Power Spectral Denstity")
xlabel('Frequency');
ylabel('Amplitude');
legend('a=0','a=0.5','a=1');

%b)
figure;
for i=1:height(PHI_psd_matrix)
    semilogy(F_axis,PHI_psd_matrix(i,:));
    hold on;
end
hold off;
grid on;
title("Power Spectral Denstity(semilogy)")
xlabel('Frequency');
ylabel('Amplitude');
legend('a=0','a=0.5','a=1');

%A3
%b)
c1=T/(10^3);
figure;
for i=1:height(PHI_psd_matrix)
    semilogy(F_axis,PHI_psd_matrix(i,:));
    hold on;
end
yline(c1);
hold off;
grid on;
title("Power Spectral Denstity(semilogy),C=T/10^3")
xlabel('Frequency');
ylabel('Amplitude');
legend('a=0','a=0.5','a=1','c1=T/10^3');

%c(
c2=T/(10^5);
figure;
for i=1:height(PHI_psd_matrix)
    semilogy(F_axis,PHI_psd_matrix(i,:));
    hold on;
end
yline(c2);
hold off;
grid on;
title("Power Spectral Denstity(semilogy),C=T/10^5")
xlabel('Frequency');
ylabel('Amplitude');
legend('a=0','a=0.5','a=1','c1=T/10^5');

%B)
k=[0:1:2*A];
TB = 10^(-2);
FsB = TB/over;
delayed_phi_matirx = [];
%B1)
%1)
for i=1:height(phit_matrix)
    figure;
     hold on;
    for j=1:length(k)-6
        %plot(t,phit_matrix(i,:),phit_matrix(i,:));
    end
    hold off;
    grid on;
    title('phi(t) and phi(t-kT) for a='+ string(a(i)));
    xlabel('Time');
    ylabel('Amplitude');
    legend('k=0','k=1','k=2');
end

%2
%{
for i=1:height(phit_matrix)
    figure;
    n=i*3;
    hold on;
    for j=length(n)-2:length(n)
       %product = phit_matrix(i,:).*delayed_phi_matirx(j,:);
       % plot(t,product);
    end
    hold off;
    grid on;
    title('Product of phi(t) with phi(t-kT) for a='+ string(a(i)));
    xlabel('Time');
    ylabel('Amplitude');

end
%}
%ores mexri stigmh 7:10







