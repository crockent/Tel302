
% EXERCISE 3
% Emmanouil-Thomas Chatzakis
% 2021030061

% RUN PART 1 AND PART 2 SEPARATELYΣ
%%
clear;
close all;
clc
%---- Part 1----%
A = 1;
hf = 4; %half duration of the pulse in symbol periods (positive integer)
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
X_I = bits_to_4_PAM(b(1:2*N),A);
X_Q = bits_to_4_PAM(b(2*N+1:4*N),A);

%%4)
[phi,t] = srrc_pulse(T,over,hf,a);
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
X_Q_mod = -2*(X_Q_t .* sin(2*pi*F0*t_conv));


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
X_psd = abs(fftshift(fft(X_mod,Nf))*Ts).^2;
Px_F_mod = X_psd/T_total;

figure;
    plot(t_conv,X_mod);
    grid on;
    title('Κυματομορφή X(t)^{mode}');
    xlabel('Χρόνος (sec)');
    ylabel('Πλάτος');  

figure;
    semilogy(F_axis,Px_F_mod,'b');
    grid on;
    title("Power Spectral Denstity of X^{mode}");
    xlabel('Frequency');
    ylabel('Amplitude');   

%8))
SNR = 20;
variance = (10*(A^2))/(Ts*(10^(SNR/10)));
gaussian_noise = sqrt(variance)*randn(1,length(X_mod));

% Adding noise 
Xmod_noise = X_mod + gaussian_noise;

%9))

X_I_demod = Xmod_noise .* cos(2*pi*F0*t_conv);
X_Q_demod = Xmod_noise .* (-sin(2*pi*F0*t_conv));

XIdemod_psd = abs(fftshift(fft(X_I_demod,Nf))*Ts).^2;
XQdemod_psd = abs(fftshift(fft(X_Q_demod,Nf))*Ts).^2;

PxI_F_demod = XIdemod_psd / T_total;
PxQ_F_demod = XQdemod_psd / T_total;


figure;
subplot(2, 1, 1);
plot(t_conv, X_I_demod);
grid on;
title('Demodulated Waveform X_I(t)');
xlabel('Time (sec)');
ylabel('Amplitude');

subplot(2, 1, 2);
plot(t_conv, X_Q_demod);
grid on;
title('Demodulated Waveform X_Q(t)');
xlabel('Time (sec)');
ylabel('Amplitude');

figure;
subplot(2, 1, 1);
semilogy(F_axis, PxI_F_demod, 'r');
grid on;
title('Power Spectral Density of Demodulated X_I');
xlabel('Frequency');
ylabel('Amplitude');

subplot(2, 1, 2);
semilogy(F_axis, PxQ_F_demod, 'k');
grid on;
title('Power Spectral Density of Demodulated X_Q');
xlabel('Frequency');
ylabel('Amplitude');

%10))
%filtering the wave forms
X_I_filtered = conv(X_I_demod,phi)*Ts;
X_Q_filtered = conv(X_Q_demod,phi)*Ts;

% Time vector for the filtered signals
t_filtered = (min(t_conv) + min(t)):Ts:(max(t_conv) + max(t));
T_total_filtered = max(t_filtered)- min(t_filtered)+1;

XI_filtered_psd = abs(fftshift(fft(X_I_filtered, Nf)) * Ts).^2;
XQ_filtered_psd = abs(fftshift(fft(X_Q_filtered, Nf)) * Ts).^2;

PxI_F_filtered = XI_filtered_psd / T_total_filtered;
PxQ_F_filtered = XQ_filtered_psd / T_total_filtered;

figure;
subplot(2, 1, 1);
plot(t_filtered, X_I_filtered);
grid on;
title('Filtered Demodulated Waveform X_I(t)');
xlabel('Time (sec)');
ylabel('Amplitude');

subplot(2, 1, 2);
plot(t_filtered, X_Q_filtered);
grid on;
title('Filtered Demodulated Waveform X_Q(t)');
xlabel('Time (sec)');
ylabel('Amplitude');

figure;
subplot(2, 1, 1);
semilogy(F_axis, PxI_F_filtered, 'r');
grid on;
title('Power Spectral Density of Filtered X_I');
xlabel('Frequency');
ylabel('Amplitude');

subplot(2, 1, 2);
semilogy(F_axis, PxQ_F_filtered, 'k');
grid on;
title('Power Spectral Density of Filtered X_Q');
xlabel('Frequency');
ylabel('Amplitude');

%11))
r=1;
for j = 2*A*over+1:over:length(t_filtered)-2*A*over
    X_I_sampled(r) = X_I_filtered(j);
    X_Q_sampled(r) = X_Q_filtered(j);
    r = r+1;
end

X_sampled = [X_I_sampled , X_Q_sampled];


scatterplot(X_sampled);
title('Scatter Plot of the Sampled Sequence');
xlabel('In-Phase Component');
ylabel('Quadrature Component');
grid on;


%12))
detected_I = detect_4_PAM(X_I_sampled, A);
detected_Q = detect_4_PAM(X_Q_sampled, A);

detected_output = [detected_I, detected_Q];


total_symbol_errors=0;
%13))
for j=1:length(X_I)
    if(X_I(j) ~= detected_I(j)||X_Q(j) ~= detected_Q(j))
    total_symbol_errors = total_symbol_errors +1;   
    end
end 

disp(['Total Symbol Errors: ', num2str(total_symbol_errors)]);

%14))
est_bits =  PAM_4_to_bits(detected_output,A,b);

%15))
bits_error_num = 0;
for j=1:length(b)
    if(b(j)~=est_bits(j))
    bits_error_num = bits_error_num +1;
    end
end
disp(['Number of bit errors: ' num2str(bits_error_num)]);

%%

%---- Part 2----%
clear;
close all;
clc;

% Parameters
A = 1;
hf = 4; % Half duration of the pulse in symbol periods (positive integer)
N = 200;
over = 10;
T = 0.01;
Ts = T / over;
a = 0.5;
Nf = 2048;
Fs = 1 / Ts;
F0 = 200;
SNRdB = 0:2:16;
K = 200; 



% Initialize error count arrays
symbol_errors = zeros(length(SNRdB), 1);
bit_errors = zeros(length(SNRdB), 1);
total_Theoretical_symbol_errors=zeros(length(SNRdB), 1);
total_Theoretical_bit_errors=zeros(length(SNRdB), 1);
for idx = 1:length(SNRdB)
    total_symbol_errors = 0;
    total_bit_errors = 0;
    SNR = 10^(SNRdB(idx) / 10); % Convert the current SNRdB value to a scalar SNR
    for r = 1:K
        % Generate random bits
        b = (sign(randn(1, 4 * N)) + 1) / 2;

        % 3) Map bits to symbols
        X_I = bits_to_4_PAM(b(1:2 * N), A);
        X_Q = bits_to_4_PAM(b(2 * N + 1:4 * N), A);

        % 4) Generate signal
        [phi, t] = srrc_pulse(T, over, hf, a);
        X_delta_I = 1 / Ts * upsample(X_I, over);
        X_delta_Q = 1 / Ts * upsample(X_Q, over);
        t_delta = [0:Ts:(N*T)-Ts];

        X_I_t = conv(X_delta_I, phi) * Ts;
        X_Q_t = conv(X_delta_Q, phi) * Ts;
        t_conv = t(1)+t_delta(1):Ts:t(end)+t_delta(end);

        % 5) Modulate the signal
        X_I_mod = 2 * (X_I_t .* cos(2 * pi * F0 * t_conv));
        X_Q_mod = -2 * (X_Q_t .* sin(2 * pi * F0 * t_conv));
        X_mod = X_I_mod + X_Q_mod;

        % 8) Add noise
        variance = (10 * (A ^ 2)) / (Ts * (10 ^ (SNR / 10)));
        gaussian_noise = sqrt(variance) * randn(1, length(X_mod));
        Xmod_noise = X_mod + gaussian_noise;

        % 9) Demodulate the noisy signal
        X_I_demod = Xmod_noise .* cos(2 * pi * F0 * t_conv);
        X_Q_demod = Xmod_noise .* sin(2 * pi * F0 * t_conv);

        % 10) Filter the demodulated signal
        X_I_filtered = conv(X_I_demod, phi) * Ts;
        X_Q_filtered = conv(X_Q_demod, phi) * Ts;
        t_filtered = (min(t_conv) + min(t)):Ts:(max(t_conv) + max(t));

        % 11) Sample the filtered signal
        r=1;
        for j = 2*A*over+1:over:length(t_filtered)-2*A*over
            X_I_sampled(r) = X_I_filtered(j);
            X_Q_sampled(r) = X_Q_filtered(j);
            r = r+1;
        end
        X_sampled = [X_I_sampled, X_Q_sampled];


        % 12) Detect the symbols
        detected_I = detect_4_PAM(X_I_sampled, A);
        detected_Q = detect_4_PAM(X_Q_sampled, A);

        % 13) Calculate symbol errors
        for j = 1:length(X_I)
            if (X_I(j) ~= detected_I(j) ||X_Q(j) ~= detected_Q(j))
                total_symbol_errors = total_symbol_errors + 1;
            end
        end

        % 14) Convert detected symbols to bits
        detected_output = [detected_I, detected_Q];
        est_bits = PAM_4_to_bits(detected_output, A, b);

        % 15) Calculate bit errors
        for j = 1:length(b)
            if b(j) ~= est_bits(j)
                total_bit_errors = total_bit_errors + 1;
            end
        end
    end

    % Experimental errors
    symbol_errors(idx) = total_symbol_errors / (K*2*length(X_I));
    bit_errors(idx) = total_bit_errors / (K*length(b));

    % Theoretical errors
    total_Theoretical_symbol_errors(idx) = (3/2)*Q(sqrt(SNR/5));
    total_Theoretical_bit_errors(idx) = total_Theoretical_symbol_errors(idx) / 4;
end


% Plot symbol error probability
figure;
semilogy(SNRdB, symbol_errors, 'r-o');
hold on;
semilogy(SNRdB, total_Theoretical_symbol_errors, 'b-s');
title('Error Probability vs. SNR');
xlabel('SNR (dB)');
ylabel('Error Probability');
legend('Symbol Error Probability', 'Theoretical Symbol Error Probability');
grid on;

%plot bit error probability
figure;
semilogy(SNRdB,bit_errors, 'r-o');
hold on;
semilogy(SNRdB, total_Theoretical_bit_errors, 'b-s');
title('Error Probability vs. SNR');
xlabel('SNR (dB)');
ylabel('Error Probability');
legend('Bit Error Probability', 'Theoretical Bit Error Probability');
grid on;















    