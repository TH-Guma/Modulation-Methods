clc; clear; close all;

N0 = 10^-2;
EbN0_dB = 0:2:6; EbN0 = 10.^(EbN0_dB /10);
Eb = EbN0 * N0; % The energy of one bit

Ntry = 5*10^3; % The number of transmitted bits
P_error_simul = zeros(1,length(EbN0_dB));
P_error_theo = zeros(1,length(EbN0_dB));

for j = 1:length(EbN0_dB)
    ts = 1/1000; % The sample time
    Tb = 1; % The time of 1 bit
    Ts = 2*Tb; % The time of 1 symbol
    V = sqrt(Eb(j)*2 /Tb); % The amplitude of the carrier signal %có thể fix Tb=Ts
    Tc = Ts/10; % The cycle of the carrier
    fc = 1/Tc; % The frequency of the carrier
    t_1symbol = 0 : ts : Ts-ts;
    s1 = V * cos(2*pi*fc*t_1symbol); % s1(t)
    s2 = V * sin(2*pi*fc*t_1symbol); % s2(t)
    s3 = -V * cos(2*pi*fc*t_1symbol); % s3(t)
    s4 = -V * sin(2*pi*fc*t_1symbol); % s4(t)
    Es = V^2*Tb; % The energy of 1 symbol
    L = length(t_1symbol); % The number of samples of 1 symbol
    Bit = randsrc(1,Ntry,[0 1]); % The transmitted bits
    s = [];
    t = [];
    for i = 1:2:Ntry
        if [Bit(i) Bit(i+1)] == [0 0]
            s = [s s1];
        elseif [Bit(i) Bit(i+1)] == [0 1]
            s = [s s2];
        elseif [Bit(i) Bit(i+1)] == [1 1]
            s = [s s3];
        else
            s = [s s4];
        end
            t_isymbol = t_1symbol + (i - 1); % Time of i-bit
            t = [t t_isymbol];
    end
% ================= The AWGN channel
    B = 1/ts; % Bandwidth of signals
    Power_noise = N0/2 * B; % The power of noise
    w = sqrt(Power_noise) * randn(1, length(s));
% ================= The received signal
    r = s +w;
% ================= Signal recovery
    phi1 = s1/sqrt(Es); % The orthonormal function 1
    phi2 = s2/sqrt(Es); % The orthonormal function 2
    h1 = flip(phi1); % The matched filter 1
    h2 = flip(phi2); % The matched filter 2

    s11 = sqrt(Es); s12 = 0;
    s21 = 0; s22 = sqrt(Es);
    s31 = -sqrt(Es); s32 = 0;
    s41 = 0; s42 = -sqrt(Es);

    Bit_rec = [];
    for i = 1:Ntry/2
        Frame = r((i-1)*L + 1 : i*L); % Construct 1 Frame with L samples of 1 symbol
        y1 = conv(Frame, h1) * ts; % r(t) passes through the matched filter 1
        r1 = y1(L);
        y2 = conv(Frame, h2) * ts; % r(t) passes through the matched filter 2
        r2 = y2(L);
        d1 = (r1 - s11)^2 + (r2 - s12)^2; % The squared distance from r to s1
        d2 = (r1 - s21)^2 + (r2 - s22)^2; % The squared distance from r to s2
        d3 = (r1 - s31)^2 + (r2 - s32)^2; % The squared distance from r to s3
        d4 = (r1 - s41)^2 + (r2 - s42)^2; % The squared distance from r to s4
        % --------- Comparator for decision
        if d1 < d2 && d1 < d3 && d1 < d4
            Bit_rec = [Bit_rec 0 0];
        elseif d2 < d1 && d2 < d3 && d2 < d4
            Bit_rec = [Bit_rec 0 1];
        elseif d3 < d1 && d3 < d2 && d3 < d4
            Bit_rec = [Bit_rec 1 1];
        else
            Bit_rec = [Bit_rec 1 0];
        end
    end
    Bit_rec;
% ================== The bit error probability
% ------------- Simulation
    [Num, rate] = biterr(Bit, Bit_rec);
    P_error_simul(j) = rate;
% ------------- Theory
    P_error_theo(j) = qfunc(sqrt(2*Eb(j)/(N0)));
end

P_error_simul
P_error_theo

figure(1)
semilogy(EbN0_dB, P_error_theo, 'r-', 'linewidth', 1.8); hold on;
semilogy(EbN0_dB, P_error_simul, 'k*', 'markersize',8);
xlabel('Eb/N0 (dB)'); ylabel('The error probability');
legend('Theory QPSK', 'Simulation')