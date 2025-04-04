function [b,prevVariable] = B_inv(cycleCount,dt,e,prevVariable,prm)
        % if cycleCount==2
        %     prevVariable.x_2_prev=ones(3,1);
        % end
        %% input
        % x_1 = prevVariable.x_1_prev;
        % x_2 = prevVariable.x_2_prev;
        % x_3 = prevVariable.x_3_prev;
        % x_debug = prevVariable.x_debug_prev;
        % x_debug2 = prevVariable.x_debug2_prev;
        % x_debug3 = prevVariable.x_debug3_prev;
        % y_w = prevVariable.y_w_prev;


        %N_LPF
        if prm.settings.isRugekuttaMethodUse == 0
            x_1 = moritaEulerMethod(@getQ_nDxdt,cycleCount,dt,[prevVariable.x_2_prev,prevVariable.x_1_prev,zeros(3,1)],prm,0);
        end
        if prm.settings.isRugekuttaMethodUse == 1
            x_1 = moritaRungekuttaMethod(@getQ_nDxdt,cycleCount,dt,[prevVariable.x_2_prev,prevVariable.x_1_prev,zeros(3,1)],prm,0);
        end



        %N
        if prm.settings.isRugekuttaMethodUse == 0 
            y_w = moritaEulerMethod(@getNDxdt,cycleCount,dt,[prevVariable.x_1_prev,prevVariable.y_w_prev,zeros(3,1)],prm,0);
        end

        if prm.settings.isRugekuttaMethodUse == 1
            y_w = moritaRungekuttaMethod(@getNDxdt,cycleCount,dt,[prevVariable.x_1_prev,prevVariable.y_w_prev,zeros(3,1)],prm,0);
        end
         
        %M-1
        x_2 = (1*prevVariable.y_w_prev + e) ./ prm.MOperatorConstPrm; 


        %Q_D
        if prm.settings.isRugekuttaMethodUse == 0
            x_3 = moritaEulerMethod(@getQ_dDxdt,cycleCount,dt,[prevVariable.x_2_prev,prevVariable.x_3_prev,zeros(3,1)],prm,0);
        end
        if prm.settings.isRugekuttaMethodUse == 1
            x_3 = moritaRungekuttaMethod(@getQ_dDxdt,cycleCount,dt,[prevVariable.x_2_prev,prevVariable.x_3_prev,zeros(3,1)],prm,0);
        end
        x_3_dot = (x_2 - x_3) ./ (prm.lowPassFilterTimePrm.D);

        %D
        a1 = prm.peltierPrm.resistance *ones(3,1);
        a2 = 2*prm.peltierPrm.seebeck * (x_3 - (prm.peltierPrm.absoluteTemperature + prm.settings.outsideTemperature));
        a3 = [prm.ma_ca13; prm.ma_ca2; prm.ma_ca13].*x_3_dot;
        for i=1:4
            a3 = a3 + ((-1)^(i+1)) *prm.A_a(:,i).*[prm.ma_ca13; prm.ma_ca2; prm.ma_ca13] .* (x_3.^i);
        end
        b1 = -a2;
        b2 = sqrt(power(a2,2) - 4*a1.*a3);
        b3 = 2*a1;
        b = [1; 0; 1].*(b1 - b2)./b3;
 

        %確かめinv       
        % y_w=b;
        % x_debug =D_inv(cycleCount,dt,...
        %                 [prevVariable.y_w_prev,prevVariable.x_debug_prev,zeros(3,1)],...
        %                 prm,0);
        % x_debug2 = prm.lowPassFilterTimePrm.D* getD_invDxdt(1,1,[b,x_debug,zeros(3,1)],prm,0) + x_debug;
        % 
        % x_debug3=x_2 - x_debug2;
        %invD→D
        %D_tilde→invD_tilde


        %次使う用
        prevVariable.x_1_prev=x_1;
        prevVariable.x_2_prev=x_2;
        prevVariable.x_3_prev=x_3;
        % prevVariable.x_debug_prev=x_debug;
        % prevVariable.x_debug2_prev=x_debug2;
        % prevVariable.x_debug3_prev=x_debug3;
        prevVariable.y_w_prev=y_w;


end

% 
% function dxdt = getD_invDxdt(cycleCount,dt,input,prm)
%         inputCurr = input(:,1);
%         y_a_prev = input(:,2);
%         y_prev = input(:,3);
% 
%         ohmHeat = -prm.peltierPrm.resistance*[1; 0; 1].*power(inputCurr,2);
%         peltierHeat =  2*prm.peltierPrm.seebeck*[1; 0; 1].*((prm.settings.outsideTemperature+prm.peltierPrm.absoluteTemperature)- y_a_prev).*inputCurr
%         dxdt = (ohmHeat + peltierHeat)./[prm.ma_ca13; prm.ma_ca2; prm.ma_ca13];
%         % 
%         for i=1:4
%             dxdt = dxdt + ((-1)^(i)) *prm.A_a(:,i).* (y_a_prev.^i);
%         end
% 
%         if prm.settings.isInterferrence
%             dxdt = dxdt + prm.interferrenceConstPrm.y_a.*[y_a_prev(2); y_a_prev(1)+y_a_prev(3);  y_a_prev(2)]./[prm.ma_ca13; prm.ma_ca2; prm.ma_ca13] + prm.interferrenceConstPrm.y_aw.*y_prev./[prm.ma_ca13; prm.ma_ca2; prm.ma_ca13];
%         end
% 
% 
% 
% end