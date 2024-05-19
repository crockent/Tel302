clear
close all
clc

%% %%%%%%%%%%%%%%%%%%%
%---------A----------%
%%%%%%%%%%%%%%%%%%%%%%

% A1

T = .01;
over = 10;
Ts = T/over;
A = 4;
a = .5;

[phi, t] = srrc_pulse(T, over, A, a);

Fs = 1/Ts;
Nf = 2048;


f = linspace(-Fs/2,Fs/2-Fs/Nf,Nf);

Phi_en = fftshift(abs(Ts*fft(phi,Nf)).^2);


figure
semilogy(f,Phi_en,'b')
grid on


% A2

N = 100;

b = (sign(randn(N, 1)) + 1)/2;

Xn = bits_to_2pam(b);

figure
stem(Xn)
grid on


t_delta = 0:Ts:(N*T-Ts);
X_delta = (1/Ts)*upsample(Xn, over);

tconv = t(1)+t_delta(1):Ts:t(end)+t_delta(end);
x = Ts*conv(X_delta, phi);

figure
plot(tconv, x)
grid on


S_x = (var(Xn)/T)*Phi_en;



% A3

X_en = fftshift(abs(Ts*fft(x,Nf)).^2);

T_total = (N + 2*A)*T;

P_x = X_en/T_total;

figure
plot(f, P_x,'b')
grid on

figure
semilogy(f, P_x,'b')
grid on

k = 50;


buffer_X = zeros(k, Nf);

for i = 1:k

    b = (sign(randn(N, 1)) + 1)/2;
    Xn = bits_to_2pam(b);

    X_delta = (1/Ts)*upsample(Xn, over);

    x_test = Ts*conv(X_delta, phi);

    X_test = fftshift(abs(Ts*fft(x_test,Nf)).^2);

    buffer_X(i,:) = X_test;

end

Sx_hat = mean(buffer_X);
    
    
figure
semilogy(f,Sx_hat,'r')
hold on
semilogy(f,S_x, 'b')
grid on


%% A4

T = .01;
over = 10;
Ts = T/over;
A = 4;
a = .5;

[phi, t] = srrc_pulse(T, over, A, a);



Fs = 1/Ts;
Nf = 2048;

Phi_en = fftshift(abs(Ts*fft(phi,Nf)).^2);

f = linspace(-Fs/2,Fs/2-Fs/Nf,Nf);

N = 100;
b0 = (sign(randn(N/2, 1)) + 1)/2;
b1 = (sign(randn(N/2, 1)) + 1)/2;

Xn = bits_to_4pam(b0,b1);

t_delta = 0:Ts:(N/2)*T-Ts;
X_delta = (1/Ts)*upsample(Xn, over);

tconv = t(1)+t_delta(1):Ts:t(end)+t_delta(end);
x = Ts*conv(X_delta, phi);

figure
plot(tconv, x)
grid on


X_en = fftshift(abs(Ts*fft(x,Nf)).^2);
T_total = (N/2 + 2*A)*T;

P_x = X_en/T_total;


figure
plot(f, P_x, 'b')
grid on

figure
semilogy(f, P_x, 'b')
grid on

S_x = (var(Xn)/T)*Phi_en;

k = 50;


buffer_X = zeros(k, Nf);

for i = 1:k
    
    b0 = (sign(randn(N/2, 1)) + 1)/2;
    b1 = (sign(randn(N/2, 1)) + 1)/2;

    Xn = bits_to_4pam(b0,b1);

    X_delta = (1/Ts)*upsample(Xn, over);

    x_test = Ts*conv(X_delta, phi);

    X_test = fftshift(abs(Ts*fft(x_test,Nf)).^2);

    buffer_X(i,:) = X_test;

end

Sx_hat = mean(buffer_X);
    
    
figure
semilogy(f,Sx_hat,'r')
hold on
semilogy(f,S_x, 'b')
grid on

%% A5

T = 2/100;
over = 20;
Ts = T/over;
A = 4;
a = .5;

[phi, t] = srrc_pulse(T, over, A, a);

Fs = 1/Ts;
Nf = 2048;


f = linspace(-Fs/2,Fs/2-Fs/Nf,Nf);

Phi_en = fftshift(abs(Ts*fft(phi,Nf)).^2);


figure
semilogy(f,Phi_en,'b')
title(sprintf("Spectral Density Energy our SRRC pulse"))
grid on

N = 100;

b = (sign(randn(N, 1)) + 1)/2;

Xn = bits_to_2pam(b);

t_delta = 0:Ts:(N*T-Ts);
X_delta = (1/Ts)*upsample(Xn, over);

tconv = t(1)+t_delta(1):Ts:t(end)+t_delta(end);
x = Ts*conv(X_delta, phi);

S_x = (var(Xn)/T)*Phi_en;

X_en = fftshift(abs(Ts*fft(x,Nf)).^2);

T_total = (N + 2*A)*T;

P_x = X_en/T_total;

figure
plot(f, P_x,'b')
grid on

figure
semilogy(f, P_x,'b')
grid on

k = 500;



buffer_X = zeros(k, Nf);

for i = 1:k

    b = (sign(randn(N, 1)) + 1)/2;
    Xn = bits_to_2pam(b);

    X_delta = (1/Ts)*upsample(Xn, over);

    x_test = Ts*conv(X_delta, phi);

    X_test = fftshift(abs(Ts*fft(x_test,Nf)).^2);

    buffer_X(i,:) = X_test;

end

Sx_hat = mean(buffer_X);
  
    
figure
semilogy(f,Sx_hat,'r')
hold on
semilogy(f,S_x, 'b')
grid on


%% %%%%%%%%%%%%%%%%%%%
%---------B----------%
%%%%%%%%%%%%%%%%%%%%%%

F0 = 300;
t = linspace(-5,5,1000);

X = randn(1,5);
PHI = (2*pi)*rand(1,5);

figure
hold on
for i = 1:5
    Y = X(i)*cos(2*pi*F0*t + PHI(i));
    plot(t,Y)
end
grid on
hold off