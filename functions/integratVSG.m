function y = integratVSG(x,fs)
%% Function to apply a simple integrator
%% Zihan Wang, Phonetics and Speech Lab, Trinity College Dublin, 2022
% Force any negative value to be positive

y=zeros(1,length(x));

Ts=1/fs;
y(1)=Ts*x(1);

% Difference equation
for n=2:length(x)
    y(n)=(Ts*x(n))+y(n-1);
    if y(n) < 0
        y(n) = 0;
    end
end
end