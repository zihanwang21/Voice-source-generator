function params = params_transformedLF23(F0,GCI,Ee,Rd)
%% Function to generate LF parameters for synthesising 
% Based on Wang, Z., & Gobl, C. (2023). A System for Generating Voice Source Signals that Implements the Transformed LF-model Parameter Control. In Interspeech 2023
%% Zihan Wang, Phonetics and Speech Lab, Trinity College Dublin, 2023 %%
%% controlling parameters: Rd, Ee and F0.
% Input: F0 in Hz
%        GCI in seconds
%        Ee in linear scale
%        Rd
% Output: a structure with all the true LF parameters for generating the LF source
%%
F0 = F0(:)';
GCI = GCI(:)';
Ee = Ee(:)';
Rd = Rd(:)';

T0 = 1./F0;
Up = 110.*Rd.*Ee./(1000.*F0);

paramlen = length(F0);
alpha = zeros(1,paramlen);
epsi = zeros(1,paramlen);
omega = zeros(1,paramlen);
Tb = zeros(1,paramlen);
Te = zeros(1,paramlen);

% Tb(end) = 0.003;
Tb(end) = 3*T0(end)*(0.048*Rd(end) - 0.01); % 3*Ta
% Tb: The Tb of the final pulse is 3*Ta or 0.003s.
%     Afterwards, the Te of the final pulse can be calculated as 
%     Te = (1+Rk)/(2RgF0), where Rg comes from the estimated omega0_g.
%     Then the Tb of the preceding pulse is: T0 - Te.
%     (Note that T0(n) = Tb(n)+Te(n+1).)

[epsi(end),alpha(end),omega(end),Te(end),~,~,~,~,Tb(end)] = Newton_IS2023_4(Ee(end),Rd(end),F0(end),Tb(end),1);
Tb(end-1) = T0(end)-Te(end);

for n = paramlen-1:-1:2
        [epsi(n),alpha(n),omega(n),Te(n)] = Newton_IS2023_4(Ee(n),Rd(n),F0(n),Tb(n));
    Tb(n-1) = T0(n)-Te(n);
end

    [epsi(1),alpha(1),omega(1),Te(1)] = Newton_IS2023_4(Ee(1),Rd(1),F0(1),Tb(1));
%% output structure
params.GCI = GCI(:);
params.Te = Te(:);
params.OQ = Te(:).*F0(:);
params.Ee = Ee(:);
params.omega = omega(:);
params.alpha = alpha(:);
params.epsi = epsi(:);
params.Tb = Tb(:);
params.F0 = F0(:);
params.Rd = Rd(:);
params.Up = Up(:); 

end