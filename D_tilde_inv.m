function [y_a, prev_variable] = D_tilde_inv(cycle_count,dt,input,prev_variable,prm,isInterferrence)
         arguments
             cycle_count
             dt
             input
             prev_variable
             prm
             isInterferrence=prm.settings.isInterferrence
         end

        %D_inv
        y_a_tilde = D_inv(cycle_count,dt,...
                        input,...
                        prm,isInterferrence);
        %y_a = y_a_tilde;


        %Q_inv
        y_a_tilde_dot = (y_a_tilde - input(:,2)) /dt;
        y_a_test = y_a_tilde + prm.lowPassFilterTimePrm.D*y_a_tilde_dot;
        y_a = y_a_test;

        prev_variable.yyy = y_a_test - y_a;






        % inputCurr = input(:,1);
        % y_a_prev = input(:,2);
        % 
        % ohmHeat = -prm.peltierPrm.resistance*[1; 0; 1].*power(inputCurr,2);
        % peltierHeat =  2*prm.peltierPrm.seebeck*[1; 0; 1].*((prm.settings.outsideTemperature+prm.peltierPrm.absoluteTemperature)- y_a_prev).*inputCurr;
        % dxdt = (ohmHeat + peltierHeat)./[prm.ma_ca13; prm.ma_ca2; prm.ma_ca13];
        % % 
        % for i=1:4
        %     dxdt = dxdt + ((-1)^(i)) *prm.A_a(:,i).* (y_a_prev.^i);
        % end
        % 
        % 
        % y_a = prm.lowPassFilterTimePrm.D.*dxdt  + y_a_tilde;


end