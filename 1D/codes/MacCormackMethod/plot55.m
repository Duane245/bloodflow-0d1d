close all

figure(1)
plot([1:1000000]*dt,captureQ1)
title('Artery 1')
xlabel('seconds')
ylabel('ml/s')

figure(2)
plot([1:1000000]*dt,captureQ2)
title('Artery 8')
xlabel('seconds')
ylabel('ml/s')

figure(3)
plot([1:1000000]*dt,captureQ3)
title('Artery 37')
xlabel('seconds')
ylabel('ml/s')

figure(4)
plot([1:1000000]*dt,captureQ4)
title('Artery 54')
xlabel('seconds')
ylabel('ml/s')

% % % % % % % % % % % % % % % % % % % % % 
P1 = beta(1)*((sqrt(captureA1))-(sqrt(A0(1))));
figure(5)
plot([1:1000000]*dt,P1/1000)
title('Artery 1')
xlabel('seconds')
ylabel('Kpa')

P2 = beta(8)*((sqrt(captureA2))-(sqrt(A0(8))));
figure(6)
plot([1:1000000]*dt,P2/1000)
title('Artery 8')
xlabel('seconds')
ylabel('Kpa')

P3 = beta(37)*((sqrt(captureA3))-(sqrt(A0(37))));
figure(7)
plot([1:1000000]*dt,P3/1000)
title('Artery 37')
xlabel('seconds')
ylabel('Kpa')

P4 = beta(54)*((sqrt(captureA4))-(sqrt(A0(54))));
figure(8)
plot([1:1000000]*dt,P4/1000)
title('Artery 54')
xlabel('seconds')
ylabel('Kpa')
