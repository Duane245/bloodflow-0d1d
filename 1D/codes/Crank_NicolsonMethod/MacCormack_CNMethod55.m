function [A,Q] = MacCormack_CNMethod55(number_of_vessel,A,Q,beta,x,L,dt,A0,Rt,rho,dx,E,a,M,MA,K,t,N,init_A,init_Q,idxOutflowArtery)

% terminal-vessel lookup built once per call (was O(28) nested loops per vessel)
isOutflow = false(number_of_vessel,1);
isOutflow(idxOutflowArtery) = true;

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% for i = 1:number_of_vessel
%     lamda1end = init_Q{i}(end)/A{i}(end) + sqrt(beta(i)/2/rho*sqrt(A{i}(end)))*100;
%     lamda20 = init_Q{i}(1)/A{i}(1) - sqrt(beta(i)/2/rho*sqrt(A{i}(1)))*100;
%     
%     Qend = interp1(x{i},init_Q{i},L(i)-lamda1end*dt);
%     Aend = interp1(x{i},A{i},L(i)-lamda1end*dt);
%     Q1 = interp1(x{i},init_Q{i},-lamda20*dt);
%     A1 = interp1(x{i},A{i},-lamda20*dt);
%     W1n1 = Qend/Aend+4*sqrt(beta(i)/2/rho*sqrt(Aend))*100;
%     W10 = 4*sqrt(beta(i)/2/rho*sqrt(A0(i)))*100;
%     W20 = -4*sqrt(beta(i)/2/rho*sqrt(A0(i)))*100;
%     W2n1 = W20-Rt(i)*(W1n1-W10);
%     Aend = (W1n1-W2n1).^4/1024*(rho/beta(i))^2*1e-8;
%     Qend = Aend*(W1n1+W2n1)/2;
%     
%     %%%%%%%%%%%%%MacCormack Method%%%%%%%%
% %     init_A = A{i};
%     % predictor step
%     Qfd = (init_Q{i}(3:N(i))-init_Q{i}(2:N(i)-1))/dx;
%     
%     dA_p = -Qfd;
%     
%     A{i}(2:N(i)-1) = A{i}(2:N(i)-1) + dA_p*dt;
%     
%     %Corrector step
%     Qrd = (init_Q{i}(2:N(i)-1)-init_Q{i}(1:N(i)-2))/dx;
%     
%     dA_c = -Qrd;
%     
%     % average time step
%     dA_avg = 0.5*(dA_p + dA_c);
%     
%     A{i}(2:N(i)-1) = init_A{i}(2:end-1) + dA_avg*dt;
%     
% %     A{i}(1) = 2*A{i}(2)-A{i}(3);
%     if i == 1
%         A{i}(1) = 2*A{i}(2)-A{i}(3);
%     end
% %     A{i}(end) = Aend;
%     
%     for j = idxOutflowArtery
%         if i == j
%             A{i}(end) = Aend;
%         end
%     end
%     
%     %%%%%%%%%%%%%%Crank_Nicolson%%%%%%%%%%%
%     K{i} = ((1-E(i))*init_Q{i}(2:end-1)+(E(i)/2)*(init_Q{i}(3:end)+init_Q{i}(1:end-2)))' + MA{i}*A{i}' + MA{i}*init_A{i}';
%     %     K(1) = K(1)-a*Q1;
% %     K{i}(1) = K{i}(1)-a(i)*InflowFuncation(t);
%     if i == 1
%         K{i}(1) = K{i}(1)-a(i)*InflowFuncation(t);
%     else
%         K{i}(1) = K{i}(1)-a(i)*Q{i}(1);
%     end
% %     K{i}(end) = K{i}(end)-a(i)*Qend;
%     K{i}(end) = K{i}(end)-a(i)*Q{i}(end);
%     
%     Q{i}(2:end-1) = M{i}\K{i};
%     
%     
%     % Boundary values - inflow
%     if i == 1
%         Q{i}(1) = InflowFuncation(t);
%     end
%     
%     for j = idxOutflowArtery
%         if i == j
%             Q{i}(end) = Qend;
%         end
%     end
%     
% %     Q{i}(end) = Qend;
%     
% 
% end
% % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
for i = 1:number_of_vessel
    % Terminal R-only boundary state: only needed for outflow vessels;
    % skipping for the 27 non-outflow vessels avoids wasted interpolation
    % and characteristic algebra.
    if isOutflow(i)
        Ai_end = A{i}(end); Qi_end = Q{i}(end);
        lamda1end = Qi_end/Ai_end + sqrt(beta(i)/2/rho*sqrt(Ai_end))*100;

        % inline 2-pt linear interp on uniform grid (replaces interp1)
        xi = x{i}; dxi = xi(2); ni = numel(xi);
        xq = L(i) - lamda1end*dt;
        pos = xq/dxi + 1; i0 = floor(pos);
        if i0 < 1, i0 = 1; elseif i0 > ni-1, i0 = ni-1; end
        fr  = pos - i0;
        Aend = (1-fr)*A{i}(i0) + fr*A{i}(i0+1);
        Qend = (1-fr)*Q{i}(i0) + fr*Q{i}(i0+1);

        ki = sqrt(beta(i)/(2*rho))*100;
        W1n1 = Qend/Aend + 4*ki*Aend^0.25;
        W10  =  4*ki*A0(i)^0.25;
        W20  = -W10;
        W2n1 = W20 - Rt(i)*(W1n1 - W10);
        Aend = (W1n1-W2n1)^4/1024*(rho/beta(i))^2*1e-8;
        Qend = Aend*(W1n1+W2n1)/2;
    end
    
    %%%%%%%%%%%%%MacCormack Method%%%%%%%%
%     init_A = A{i};
    % predictor step
    Qfd = (Q{i}(3:N(i))-Q{i}(2:N(i)-1))/dx;
    
    dA_p = -Qfd;
    
    A{i}(2:N(i)-1) = A{i}(2:N(i)-1) + dA_p*dt;
    
    %Corrector step
    Qrd = (Q{i}(2:N(i)-1)-Q{i}(1:N(i)-2))/dx;
    
    dA_c = -Qrd;
    
    % average time step
    dA_avg = 0.5*(dA_p + dA_c);
    
    A{i}(2:N(i)-1) = init_A{i}(2:N(i)-1) + dA_avg*dt;
    
%     A{i}(1) = 2*A{i}(2)-A{i}(3);
    if i == 1
        A{i}(1) = 2*A{i}(2)-A{i}(3);
    end
%     A{i}(end) = Aend;
    if isOutflow(i)
        A{i}(end) = Aend;
    end
    
    %%%%%%%%%%%%%%Crank_Nicolson%%%%%%%%%%%
    K{i} = ((1-E(i))*init_Q{i}(2:end-1)+(E(i)/2)*(init_Q{i}(3:end)+init_Q{i}(1:end-2)))' + MA{i}*A{i}' + MA{i}*init_A{i}';
    %     K(1) = K(1)-a*Q1;
%     K{i}(1) = K{i}(1)-a(i)*InflowFuncation(t);
    if i == 1
        K{i}(1) = K{i}(1)-a(i)*InflowFuncation(t);
        
    else
        K{i}(1) = K{i}(1)-a(i)*Q{i}(1);
    end
%     K{i}(end) = K{i}(end)-a(i)*Qend;
%     K{i}(end) = K{i}(end)-a(i)*Q{i}(end);
    
    if isOutflow(i)
        K{i}(end) = K{i}(end)-a(i)*Qend;
    else
        K{i}(end) = K{i}(end)-a(i)*Q{i}(end);
    end
    
    Q{i}(2:end-1) = M{i}\K{i};
    
    
    % Boundary values - inflow
    if i == 1
        Q{i}(1) = InflowFuncation(t);
    end
    
    if isOutflow(i)
        Q{i}(end) = Qend;
    end
    
%     Q{i}(end) = Qend;
    

end
end
