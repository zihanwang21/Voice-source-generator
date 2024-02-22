function S_ah = aspNoiseVSG(pulse,F0,Ee,Up,AH_dB,ratio,fs)
%% Generate aspiration noise based on Gobl (2006).
%% Zihan Wang, Phonetics and Speech Lab, Trinity College Dublin, 2022
% for a single pulse only

% coding notes:
% AH = .000000023;
AH = db2mag(AH_dB);

pulse_int = integratVSG(pulse,fs);
% Note that integrat1 function has a built-in compensation for fs 

% Up = (110*Rd*Ee)/(1000*F0);
% Rd = 1000*Up*F0/(110*Ee);
% T0 = 1/F0;

% U_ac = (42*(Ee/Up)*T0)-91;
% Theoretically this is supposed to be set and used for calculating the ratio.
U_dc = (83*(Up/Ee))+34;

S_ah = (Ee^1.35).*(F0^1.05).*sqrt((ratio.*pulse_int+U_dc)).*AH.*randn(1,length(pulse));
end

