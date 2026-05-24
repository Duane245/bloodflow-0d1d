function Qin = InflowFuncation(t)
% Qc = 1;
Qc = 1;
Tc = 0.4;
% Qin = Qc*sin((2*pi/Tc)*t)*HeavisideFuncation(Tc/2-t);

t0 = 0.02;
tau = 0.03;
Qin = Qc*exp(-4*pi*(t-t0)^2/tau^2);
end
