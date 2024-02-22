function [epsi,alpha,omega,Te,cnt,f1,f2,kk,Tb] = Newton_IS2023_4(Ee,Rd,F0,Tb,indicator_TbUpdate)
%% Newton-Raphson Method for calculating epsi0lon, alpha and Te %%
% For more info see Wang, Z. and Gobl, C. A System for Generating Voice
% Source Signals that Implements the Transformed LF-model Parameter
% Control. Proc. INTERSPEECH 2023, Dublin, Ireland, 2023.
%% Zihan Wang, Phonetics and Speech Lab, Trinity College Dublin, 2023
% NOTE: For one pulse only. 
% Ee in linear scale; F0 as Hz; Tb in seconds.
% UPDATE Tb every iteration
%% extended Rd range with kk
% Rd~(0,1.71]: kk = 1;
% Rd~(1.71,2.7]: kk = 0.031978*Rd^3-0.32428*Rd^2+1.0706*Rd-0.044299;
% Rd~(1.71,2.7]: kk = ;
% Rd~(1.71,2.7]: kk = ;
%% initial settings %%
iter_max = 999;
ep = 1e-10; % error precision
cnt = 1;

if nargin<5
    indicator_TbUpdate = 0;
end

Up = 110*Rd*Ee/(1000*F0);
T0 = 1/F0;

change_epsi = 1;
change_alpha = 1;
change_Te = 0.1*T0;
%% Extended predictions of Ra & Rk - Huber et al. 2012, 2013 %%
if Rd<=2.7
    Rk = 0.118*Rd + 0.224;
    if Rd<0.21
        Ra = 1e-3; 
    else
        Ra = 0.048*Rd - 0.01;
    end
else
    %% Fant's approach:
    Rk = 1 - 1.230940525/(Rd+2/(2.17*5.96)); 
    Ra = 0.32292/Rd; 
end

kk = 1;
if Rd >1.71 && Rd<=2.7
    kk = 0.0319776606698207*Rd^3-0.324280056988953*Rd^2+1.07061117027544*Rd-0.0442987850160119;
elseif Rd>2.7 && Rd<=3.7
    kk = 0.481*(log(Rd))^3 - 1.9389*(log(Rd))^2 + 2.5859*log(Rd) - 0.0664;
elseif Rd>3.7 && Rd<=5.7
    kk = exp(1/(-0.14228292*Rd^3 + 2.3586596*Rd^2 -10.36752*Rd + 27.053424));
elseif Rd>5.7
    kk = exp(1/(3.2274*Rd - 0.155));
end
% kk = kk+0.01;
Rk = Rk*kk;

Ta = T0*Ra;

if Tb<=Ta
    Tb = Ta*1.001;
end

%% initial vals
epsi0 = 1./Ta;

% if Rd>0.41
%     Te0 = 0.5*T0;
% else
%     Te0 = Rd*T0;
% end

if Rd<=1.71
    Te0 = 0.3*Rd*T0;
    alpha0 = 1;
%     indicator_TbUpdate = 1;
%     Te0_max = T0 - Ta;
else
    indicator_TbUpdate = 0; % Tb will not be updated in all cases when Rd>2.7
    if Rd<=2.7
        b2 = -0.11669148; b1 = -0.08011759; b0 = 2.178022204;
        c2 = 0.095942248; c1 = 0.613891367; c0 = -4.133302124;
        OQ_min = -(c2*Rd^2+c1*Rd+c0)/(2*(b2*Rd^2+b1*Rd+b0));
        alpha0 = exp(((-42.1*(log(OQ_min))^2 +10.4*log(OQ_min) - 3.79)*Ra + 2.61*log(OQ_min) - 2.17)*Rd +  exp((-6.75*log(OQ_min) - 0.55)*Ra + 2.16));
    else
        OQ_min_27 = 0.799432204514675;
        OQ_min = OQ_min_27 + (1 - OQ_min_27)/(0.32292/2.7)*(0.32292/2.7 - 0.32292/Rd);
        alpha0 = 0;
    end
    Te0 = OQ_min*T0;
%     Te0_max = OQ_min*T0;
    Tb = (1-OQ_min)*T0;
end

% Tb = T0-Te0;
alpha_p = alpha0;

alpha0_max = 50000;
alpha0_min = -10000;
Te0_max = T0 - Ta;
Te0_min = 0.01*T0;

%% iterations %%
while cnt<=iter_max && (abs(change_epsi)>ep||(abs(change_alpha)>ep || abs(change_Te)>ep))
    % epsi0
    f_epsi0 = epsi0.*Ta - 1 + exp(-epsi0.*Tb);
    f_epsi0_der = Ta - Tb.*exp(-epsi0.*Tb);
    change_epsi = f_epsi0./f_epsi0_der;
    epsi0 = epsi0 - change_epsi;

    if epsi0 < 0.0001
        epsi0 = 0.0001; % else when epsilon gets too small the frequency domain LF does not work
    end
    
    S1 = pi*(1+Rk);
    S2 = exp(-alpha0*Te0);
    S3 = exp(-alpha0*Te0*Rk/(1+Rk));
    S4 = 1 - exp(-epsi0*Tb);
    S5 = 1 + alpha0*Te0;
    S6 = 1 - exp(-epsi0*Tb)*(1+epsi0*Tb);
    
    f1 = Up*((S1/Te0)^2+alpha0^2)*sin(S1) + Ee*S1/Te0*(S2+S3);
    f2 = epsi0*S4*(alpha0*sin(S1)-S1*cos(S1)/Te0+S1*S2/Te0)+S6*((S1/Te0)^2+alpha0^2)*sin(S1);

    %% partial derivatives %%
    par_f1_Te = -2*S1^2*Te0^-3*Up*sin(S1) - Ee*pi*Te0^-2*(S5*(1+Rk)*S2+(1+Rk*S5)*S3);
    par_f1_alpha = 2*alpha0*Up*sin(S1) - pi*Ee*((1+Rk)*S2+Rk*S3);
    par_f2_Te = epsi0*S4*S1*Te0^-2*(cos(S1)-S5*S2) - 2*S1^2*Te0^-3*sin(S1)*S6;
    par_f2_alpha = epsi0*S4*(sin(S1)-S1*S2) + 2*alpha0*sin(S1)*S6;

    det_J = par_f1_Te*par_f2_alpha - par_f1_alpha*par_f2_Te;
    numerator_Te = par_f2_alpha*f1 - par_f1_alpha*f2;
    numerator_alpha = par_f1_Te*f2 - par_f2_Te*f1;

    change_Te = numerator_Te/det_J;
    change_alpha = numerator_alpha/det_J;

    alpha0 = alpha0 - change_alpha;
    Te0 = Te0 - change_Te;
    
    cnt = cnt+1;

    %% constraints %%
    
    if Te0 > Te0_max
        Te0 = Te0_max;
    elseif Te0 < 0
        Te0 = - Te0;
    end

    if alpha0 > alpha0_max
        alpha0 = alpha0_max - mod(alpha0 - alpha0_max,alpha0_max - alpha_p);
    elseif alpha0 < alpha0_min
        alpha0 = alpha0_max - mod(alpha0_min - alpha0,alpha0_max - alpha_p);
    end

    if indicator_TbUpdate==1
        Tb = T0-Te0;
    end
end

omega0 = pi*(1+Rk)/Te0;

if cnt==iter_max+1
    warning(append('Did not converge: Rd = ',num2str(Rd),' Ee = ',num2str(Ee),' F0 = ',num2str(F0),' Tb = ',num2str(Tb)))
%     Te0 = OQ_min*T0;
%     [epsi0,alpha0,omega0,cnt] = Newton_Gobl2017_3(Ee,Rd,F0,Te0,Tb);
end

epsi = epsi0;
alpha = alpha0;
Te = Te0;
omega = omega0;

end