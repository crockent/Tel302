clear;
close all;
clc
% EXERCISE 3
% Emmanouil-Thomas Chatzakis
% 2021030061

%%Part 1%%
A = 1;
N=200;
over =10;
T = 0.01;
Ts = T/over;
a= 0.5;
Nf = 2048;
Fs = 1/Ts;
F0 = 200;
%%1)
b= (sign(randn(1,4*N))+1)/2;

%%3)
X_I = bits_to_4_PAM(b(1:2*N),A)
X_Q = bits_to_4_PAM(b(2*N+1:4*N),A)

%%4)
[phi,t] = srrc_pulse(T,over,A,a);
X_delta_I = 1/Ts*upsample(X_I, over);
X_delta_Q = 1/Ts*upsample(X_Q, over);

t_delta = [0:Ts:(N*T)-Ts];


X_I_t =  conv(X_delta_I,phi)*Ts;
X_Q_t =  conv(X_delta_Q,phi)*Ts;
t_conv = t(1)+t_delta(1):Ts:t(end)+t_delta(end);

figure;
    subplot(2, 1, 1);
    plot(t_conv,X_I_t);
    grid on;
    title('Κυματομορφή X_I(t)');
    xlabel('Χρόνος (sec)');
    ylabel('Πλάτος');
    
    
    subplot(2, 1, 2);
    plot(t_conv, X_Q_t);
    grid on;
    title('Κυματομορφή X_Q(t)');
    xlabel('Χρόνος (sec)');
    ylabel('Πλάτος');

   
XI_psd = abs(fftshift(fft(X_I_t,Nf))*Ts).^2;
XQ_psd = abs(fftshift(fft(X_Q_t,Nf))*Ts).^2;
F_axis = Fs*linspace(-1/2,1/2-Fs/Nf,Nf);

T_total = max(t_conv)- min(t_conv)+1;

PxI_F = XI_psd/T_total;
PxQ_F = XQ_psd/T_total;

figure;
    subplot(2,1,1);
    semilogy(F_axis,PxI_F,'r');
    grid on;
    title("Power Spectral Denstity of X_I");
    xlabel('Frequency');
    ylabel('Amplitude');

    subplot(2,1,2);
    semilogy(F_axis,PxQ_F,'k');
    grid on;
    title("Power Spectral Denstity of X_Q");
    xlabel('Frequency');
    ylabel('Amplitude');


%%5)
X_I_mod = 2*(X_I_t .* cos(2*pi*F0*t_conv));
X_Q_mod = -2*(X_Q_t .* cos(2*pi*F0*t_conv));


figure;
    subplot(2, 1, 1);
    plot(t_conv,X_I_mod);
    grid on;
    title('Κυματομορφή X_I(t)^{mode}');
    xlabel('Χρόνος (sec)');
    ylabel('Πλάτος');
    
    
    subplot(2, 1, 2);
    plot(t_conv, X_Q_mod);
    grid on;
    title('Κυματομορφή X_Q(t)^{mode}');
    xlabel('Χρόνος (sec)');
    ylabel('Πλάτος');

XImode_psd = abs(fftshift(fft(X_I_mod,Nf))*Ts).^2;
XQmode_psd = abs(fftshift(fft(X_Q_mod,Nf))*Ts).^2;




PxI_F_mod = XImode_psd /T_total;
PxQ_F_mod = XQmode_psd/T_total;

figure;
    subplot(2,1,1);
    semilogy(F_axis,PxI_F_mod,'r');
    grid on;
    title("Power Spectral Denstity of X_I^{mode}");
    xlabel('Frequency');
    ylabel('Amplitude');

    subplot(2,1,2);
    semilogy(F_axis,PxQ_F_mod,'k');
    grid on;
    title("Power Spectral Denstity of X_Q^{mode}");
    xlabel('Frequency');
    ylabel('Amplitude');

  %%6)
X_mod = X_I_mod + X_Q_mod;
figure;
    plot(t_conv,X_mod);
    grid on;
    title('Κυματομορφή X(t) mode');
    xlabel('Χρόνος (sec)');
    ylabel('Πλάτος');  

X_psd = abs(fftshift(fft(X_mod,Nf))*Ts).^2;
Px_F_mod = X_psd/T_total;

figure;
    semilogy(F_axis,Px_F_mod,'b');
    grid on;
    title("Power Spectral Denstity of X^{mode}");
    xlabel('Frequency');
    ylabel('Amplitude');   

















    