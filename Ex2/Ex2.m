%Emmanouil-Thomas Chatzakis 2021030061
clear all; close all; clc
%A1
T= 10^(-3);
over = 10;
Ts = T/over;
A=4;
a=0.5;

[phi,t] = srrc_pulse(T,over,A,a);

Nf=2048;
Fs = 1/Ts;
f_axis = (-0.5:1/Nf:0.5-1/Nf);
F_axis = Fs*f_axis;

PHI = fftshift(fft(phi,Nf))*Ts;
XF_abs = abs(PHI);
PHI_psd = XF_abs.^2;

figure;
semilogy(F_axis,PHI_psd,'k');
grid on;
title("Power Spectral Denstity");
xlabel('Frequency');
ylabel('Amplitude');

%A2
N=100;

b= (sign(randn(1,N))+1)/2;
X = bits_to_2pam(b);

figure;
stem(X);
grid on;

X_delta = 1/Ts*upsample(X, over);
t_delta = [0:Ts:(N*T)-Ts];

X_t = conv(X_delta,phi)*Ts;
t_conv = [min(t_delta)+min(t):Ts:max(t_delta)+max(t)];

figure;
plot(t_conv,X_t,'r');
grid on;

S_x = (var(X)/T)*PHI_psd;

%A3
X_F = fftshift(fft(X_t,Nf))*Ts;
XF_abs = abs(X_F);
XF_psd = XF_abs.^2;

%T_total = (N + 2*A)*T;
T_total = max(t_conv)- min(t_conv);

Px_F = XF_psd/T_total;

figure;
plot(F_axis,Px_F,'r');
grid on;

figure;
semilogy(F_axis,Px_F,'r');
grid on;

k=500;

X_tests = zeros(k,Nf);
for i=1:k
    b_test= (sign(randn(1,N))+1)/2;
    X_test = bits_to_2pam(b_test);
    X_delta_test = 1/Ts*upsample(X_test, over);
    X_t_test = conv(X_delta_test,phi)*Ts;
    X_F_test = fftshift(fft(X_t_test,Nf))*Ts;
    XF_abs_test = abs(X_F_test);
    XF_psd_test = XF_abs_test.^2;
    X_tests(i,:)=XF_psd_test;
end

Sx_tests = mean(X_tests);

figure;
semilogy(F_axis,S_x,'b');
hold on;
semilogy(F_axis,Sx_tests,'r');
hold off;

%A4

b1= (sign(randn(1,N/2))+1)/2;
b2= (sign(randn(1,N/2))+1)/2;

X4 = bits_to_4pam(b1,b2);

X4_delta=1/Ts*upsample(X4, over);
t_delta4 = [0:Ts:(N/2*T)-Ts];

X4_t = conv(X4_delta,phi)*Ts;
t_conv4 = [min(t_delta4)+min(t):Ts:max(t_delta4)+max(t)];

figure;
plot(t_conv4,X4_t);
grid on;

S_x4 = (var(X4)/T)*PHI_psd;

%A5




%Β
clear all; close all; clc
A=4;
F0 = 500;
tB = linspace(-A,A,2*F0 );
XB= randn(1,5);
fi = 2*pi*rand(1,5);

figure;
hold on;
for i=1:5
Y = XB(i)*cos(2*pi*F0*tB+fi(i));
plot(tB,Y);

end
grid on;
hold off;
xlabel("time axis");
ylabel("Amplitude");
title("Realizations of Y");

