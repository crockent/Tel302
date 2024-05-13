clear all; close all; clc
%A1
T=10^-2;
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
delayed_phi_matrix = [];
%B1)
%a
for i=1:height(phit_matrix)
    figure;
     hold on;
    for j=1:length(k)-6
       delayed_phi = delayseq(phit_matrix(i,:).',k(j)*T,Fs).';
       delayed_phi_matrix = [delayed_phi_matrix; delayed_phi];
       plot(t,delayed_phi);
    end
    hold off;
    grid on;
    title('phi(t) and phi(t-kT) for a='+ string(a(i)));
    xlabel('Time');
    ylabel('Amplitude');
    legend('k=0','k=1','k=2');
end

%b
for i=1:height(phit_matrix)
    for j=1:length(k)-6
        product = phit_matrix(i,:).*delayed_phi_matrix(j,:);
        figure;
        plot(t,product); 
        title("phi(t)*phi(t-kT) for" + " " + "k=" + string(k(j))+ " " + "and a=" + " " +string(a(i)));
        xlabel('Time');
        ylabel('Amplitude');
        grid on;
    end
end



for i=1:height(phit_matrix)
    for j=1:length(k)-5
       delayed_phi_2 = delayseq(phit_matrix(i,:).',k(j)*T,Fs).';
       integral = sum(phit_matrix(i,:).*delayed_phi_2)*Ts;
       fprintf("Integral phi(t)*phi(t-kT) for" + " " + "k=" + string(k(j))+ " " + "and a=" + " " +string(a(i)) +" " +"is:"+ " " + string( integral) +'\n');

    end
  
end

 




%C1
N1=100;


b= (sign(randn(1,N1))+1)/2;

%C2
%a
X = bits_to_2pam(b);
X_delta = 1/Ts*upsample(X, over);

%b
tN = [0:Ts:(N1*T)-Ts];
figure();
stem(tN,X_delta);
title('Delta pulse train');
xlabel('Time');
ylabel('Amplitude');
grid on;


%g
[phi_c2,t_c2] = srrc_pulse(T,over,A,a(2));
X_t = conv(X_delta,phi_c2)*Ts;


t_conv = [min(tN)+min(t_c2):Ts:max(tN)+max(t_c2)];
figure;
plot(t_conv,X_t);
title('Convolution of Delta pulse train with phi');
xlabel('Time');
ylabel('Amplitude');
grid on;

%d
phi_flipped = phi_c2(numel(phi_c2):-1:1);
t_flipped = t_c2(numel(t_c2):-1:1);
Z_t = conv(X_t,phi_flipped)*Ts;
t_conv2 = [min(t_conv)+min(t_flipped):Ts:max(t_conv)+max(t_flipped)];
figure;
plot(t_conv2,Z_t);
hold on;
stem([0:N1-1]*T,X);
hold off;
title('Z(t)');
xlabel('Time');
ylabel('Amplitude');
grid on;

p= 0:Ts
%ores mexri stigmh 13:00
fprintf('\n');
for i=1:N1
    
    fprintf("Values of symbols Vecotr X(k) for k ="+" "+string(i-1)+" "+ ":" +string(X(i))+" "+"Values of Z(kT) for k= "+" "+string(i-1)+" "+ ":" +string(Z_t(i))  +'\n');
end


