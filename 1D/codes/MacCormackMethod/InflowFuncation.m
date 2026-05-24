function Qin = InflowFuncation(t)
% Qc = 1;
Qc = 1;
Tc = 0.4;
% Tc = 0.8;
Qin = Qc*sin((2*pi/Tc)*t)*HeavisideFuncation(Tc/2-t);             % 单周期
% Qin = Qc*sin((2*pi/Tc)*t)*HeavisideFuncation(sin((2*pi/Tc)*t));     % 多周期
end
