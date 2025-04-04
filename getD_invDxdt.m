function dxdt = getD_invDxdt(cycleCount,dt,input,prm,isInterferrence)
         arguments
             cycleCount
             dt
             input
             prm
             isInterferrence=prm.settings.isInterferrence
         end

         
        inputCurr = input(:,1);
        y_a_prev = input(:,2);
        y_prev = input(:,3);

        ohmHeat = -prm.peltierPrm.resistance*[1; 0; 1].*power(inputCurr,2);
        peltierHeat =  2*prm.peltierPrm.seebeck*[1; 0; 1].*((prm.settings.outsideTemperature+prm.peltierPrm.absoluteTemperature)- y_a_prev).*inputCurr;
        dxdt = (ohmHeat + peltierHeat)./[prm.ma_ca13; prm.ma_ca2; prm.ma_ca13];
        % 
        for i=1:4
            dxdt = dxdt + ((-1)^(i)) *prm.A_a(:,i).* (y_a_prev.^i);
        end

        if isInterferrence
            dxdt = dxdt + prm.interferrenceConstPrm.y_a.*([y_a_prev(2); -y_a_prev(2)+y_a_prev(1)+y_a_prev(3);  y_a_prev(2)] - y_a_prev)./[prm.ma_ca13; prm.ma_ca2; prm.ma_ca13] ...
                        + prm.interferrenceConstPrm.y_aw.*y_prev./[prm.ma_ca13; prm.ma_ca2; prm.ma_ca13];
        end


end