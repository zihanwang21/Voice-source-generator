function [pulse,A_LF,Ph_LF,H_LF] = freq_LF_pulse(Ee,Te,omega,alpha,epsi,Tb,fs,l)
%% Function that generates one single LF pulse in frequency domain with the LF parameters.
% Based on Gobl, C. (2021). The LF Model in the Frequency Domain for Glottal Airflow Modelling Without Aliasing Distortion. In Interspeech (pp. 2651-2655).
%% Zihan Wang, Phonetics and Speech Lab, Trinity College Dublin, 2022
%%
% l = fs.*(Te+Tb);
% l = 256;
% f=(1:l)'*fs/l;
f=(0:l-1)'*fs/l;


% H_open = Ee.*(2.*pi.*f.*1i - alpha + omega.*csc(omega.*Te).*(cos(omega.*Te) - exp((2.*pi.*f.*1i-alpha).*Te)))./((2.*pi.*f.*1i - (alpha+omega)).*(2.*pi.*f.*1i - (alpha-omega)));
% H_ret = -Ee.*(2.*pi.*f.*1i - epsi.*exp(-epsi.*Tb)./(1-exp(-epsi.*Tb)).*(1-exp(-2.*pi.*f.*1i.*Tb)))./(2.*pi.*f.*1i.*(2.*pi.*f.*1i + epsi));
% H_LF = H_open+H_ret;

Re1 = omega.*(cot(omega.*Te) - exp(-alpha.*Te).*cos(2.*pi.*f.*Te)./sin(omega.*Te)) - alpha;
Im1 = 2.*pi.*f - omega.*exp(-alpha.*Te).*sin(2.*pi.*f.*Te)./sin(omega.*Te);

Re2 = epsi.*exp(-epsi.*Tb).*(1-cos(2.*pi.*f.*Tb))./(1-exp(-epsi.*Tb));
Im2 = epsi.*exp(-epsi.*Tb).*sin(2.*pi.*f.*Tb)/(1-exp(-epsi.*Tb)) - 2.*pi.*f;

Ao = Ee.*sqrt(Re1.^2 +Im1.^2)./(sqrt(alpha.^2+(2.*pi.*f - omega).^2).*sqrt(alpha.^2+(2.*pi.*f + omega).^2));
Pho = atan2(Im1,Re1) - atan2(2.*pi.*f - omega, -alpha) - atan2(2.*pi.*f + omega, -alpha);
Ar = Ee.*sqrt(Re2.^2 +Im2.^2)./(2.*pi.*f.*sqrt(epsi.^2+(2.*pi.*f).^2));
Phr = atan2(Im2,Re2) - atan2(2.*pi.*f,epsi) - pi./2;

A_LF = sqrt(Ao.^2+Ar.^2+2.*Ao.*Ar.*cos(Pho-Phr));
Ph_LF = atan2(Ao.*sin(Pho)+Ar.*sin(Phr),Ao.*cos(Pho)+Ar.*cos(Phr));
Ph_LF = Ph_LF - 2.*pi.*f.*Te; 
% Phase shift to make the pulse start at To = 0 rather than at Te


Re_LF = A_LF.*cos(Ph_LF);
Re_LF(1) = 0;
% Re_LF = [Re_LF(1:l/2);flip(Re_LF(1:l/2))];
Im_LF = A_LF.*sin(Ph_LF);
Im_LF(1) = 0;
% Im_LF = [Im_LF(1:l/2);-flip(Im_LF(1:l/2))];

H_LF = Re_LF + Im_LF.*1i;
pulse_temp = ifft(H_LF,'symmetric');
% pulse = [real(pulse_temp(l-Te*fs:l));real(pulse_temp(1:Tb*fs-1))];
% pulse_cl = pulse_temp(1:round(Tb*fs)-1);
% pulse_op = pulse_temp(l-round(Te*fs):l);
pulse = pulse_temp(1:round((Te+Tb).*fs));
end

