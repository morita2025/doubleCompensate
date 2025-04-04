function y_a = D_inv(cycleCount,dt,input,prm,isInterferrence)
         arguments
             cycleCount
             dt
             input
             prm
             isInterferrence=prm.settings.isInterferrence
         end

        %%input
        % inputCurr = input(:,1);
        % y_a_prev = input(:,2);
        % y_prev = input(:,3)


        if prm.settings.isRugekuttaMethodUse == 0
            y_a = moritaEulerMethod(@getD_invDxdt,cycleCount,dt,input,prm,isInterferrence);
        end

        if prm.settings.isRugekuttaMethodUse == 1
            y_a = moritaRungekuttaMethod(@getD_invDxdt,cycleCount,dt,input,prm,isInterferrence);
        end


end


