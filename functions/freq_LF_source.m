function [sig,H_LF] = freq_LF_source(params,siglength,fs,aspInd,AH_dB)
%% Function that generates the LF source signal in frequency domain with the given LF parameters.
% Based on Gobl, C. (2021). The LF Model in the Frequency Domain for Glottal Airflow Modelling Without Aliasing Distortion. In Interspeech 2021 (pp. 2651-2655).
%% Zihan Wang, Phonetics and Speech Lab, Trinity College Dublin, 2022 %%
% Input:
% params: structure that contains LF parameters 
%   The following values must be given:
%           - Ee (in linear scale)
%           - Te (in second)
%           - Tb (in second)
%           - GCI: glottal closure instants (in second)
%           - omega: for omega_g
%           - alpha
%           - epsi: for epsilon
%           - F0 (in Hz)
%           - Up
% siglength: the length of the original speech, aka the LF source signal
%            (in samples, i.e. duration in second * fs)
% fs: sampling frequency (in Hz)
% aspInd: indicator for aspiration noise (0: no aspNoise; 1: with
%         aspNoise), modelled following Gobl (2006)

% Output: sig: LF source signal

%% presets
if nargin < 5 
    AH_dB = -120; % default setting: -120 dB 
    if nargin < 4
        aspInd = 0; % default setting: no aspiration noise
    end
end

l = 2^(nextpow2(siglength)); % ifft length
concInd = 1; 

% concatenate method: 0 = time-domain concatenate
%                     1 = frequency-domain phasor addition

%  The concatenation described in Gobl (2021) should be realised in frequency domain. 

%  Here, the concatenation can be realised in two ways:
%  in time domain (concid = 0) - the pulses are transformmed from frequency
%  domain using IFFT and then concatenated in time domain;
%  OR in frequency domain (concid = 1) as described in Gobl (2021).

F0min=20;
F0max=500;

%% unpack parameters
Ee = params.Ee;
Te = params.Te;
omega = params.omega;
alpha = params.alpha;
epsi = params.epsi;
Tb = params.Tb;
GCI = params.GCI;
F0 = params.F0;
Up = params.Up;

aspNoise = randn(1,siglength);
[maxUp,max_ind] = max(Up);
U_ac = (42*(Ee(max_ind)/maxUp)./F0(max_ind))-91;
ratio = U_ac/maxUp;
%% concatenate

if concInd == 0
    sig = zeros(1,siglength)';
    

    for n=1:length(GCI)
        if F0(n) > F0min && F0(n) < F0max
            pulse = freq_LF_pulse(Ee(n),Te(n),omega(n),alpha(n),epsi(n),Tb(n),fs,l);
            pulse = pulse.*fs;
            aspNoise_pulse = aspNoiseVSG(pulse,F0(n),Ee(n),Up(n),AH_dB,ratio,fs);
            
%             GCI_cur = round((GCI(n).*fs));
%             start = GCI_cur - round(Te(n).*fs);
            start = round((GCI(n)-Te(n)).*fs);
            stop = start + length(pulse)-1; % avoid rounding error
%             stop = GCI_cur + round(Tb(n).*fs-1);
            if start > 0 && stop <=siglength
                sig(start:stop)=sig(start:stop)+pulse;
                aspNoise(start:stop) = aspNoise_pulse;
            end
        end
    end
    sig = sig';

    if aspInd == 1
        sig = sig + aspNoise;
    end
    
else
    start_fir = round((GCI(1)-Te(1)).*fs)-1; % t0 of the first pulse
    sig_start = zeros(start_fir,1); % unvoiced region, onset
    stop_fin = round((GCI(end)+Tb(end)).*fs); % tc of the final pulse
    sig_stop = zeros(siglength-stop_fin+1,1); % unvoiced region, end
    
    shift = 0;
    H_LF = zeros(l,1);
    
    for n=1:length(GCI)
        if F0(n) > F0min && F0(n) < F0max
            [~,A_LF,Ph_LF] = freq_LF_pulse(Ee(n),Te(n),omega(n),alpha(n),epsi(n),Tb(n),fs,l);
            f=(0:l-1)'*fs/l;
            Ph_LF = Ph_LF - 2.*pi.*f.*shift;
            
            Re_cur = A_LF.*cos(Ph_LF);
            Re_cur(1) = 0;
            Im_cur = A_LF.*sin(Ph_LF);
            Im_cur(1) = 0;
            
            H_cur = Re_cur + Im_cur.*1i;
            
            GCI_cur = round((GCI(n).*fs));
            start = GCI_cur - round(Te(n).*fs);
            stop = GCI_cur + round(Tb(n).*fs-1);
            if start > 0 && stop <=siglength
                H_LF = H_LF+H_cur;
            end
        end
        shift = shift+Te(n)+Tb(n);

    end
    sig_mid = ifft(H_LF,'symmetric');
    sig_mid = sig_mid(1:stop_fin-start_fir-1);
    sig = [sig_start;sig_mid;sig_stop];
    sig = fs.*sig';

%     Ee_actual = params.Ee;
%     for n=1:length(GCI)
%         start = round((GCI(n)-Te(n)).*fs);
%         stop = round((GCI(n)+Tb(n)).*fs-1);
%         pulse = sig(start:stop);
%         Ee_actual(n) = abs(min(pulse));
%     end

    if aspInd == 1
        
        for n=1:length(GCI)
            start = round((GCI(n)-Te(n)).*fs);
            stop = round((GCI(n)+Tb(n)).*fs-1);
            pulse = sig(start:stop);
    
            aspNoise_pulse = aspNoiseVSG(pulse,F0(n),Ee(n),Up(n),AH_dB,ratio,fs);
            aspNoise(start:stop) = aspNoise_pulse;
        end
        sig = sig + aspNoise;
    end

end


end

