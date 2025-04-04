function [y, prev_variable] = N_tilde(cycle_count,dt,input,prev_variable,prm,isInterferrence)
         arguments
             cycle_count
             dt
             input
             prev_variable
             prm
             isInterferrence=prm.settings.isInterferrence
         end

         y_a = input(:,1);
         y_a_tilde = prev_variable.y_a_tilde;

        %N_LPF
        if prm.settings.isRugekuttaMethodUse == 0
            y_a_tilde = moritaEulerMethod(@getQ_nDxdt,cycle_count,dt,[y_a,y_a_tilde,zeros(3,1)],prm,isInterferrence);
        end
        if prm.settings.isRugekuttaMethodUse == 1
            y_a_tilde = moritaRungekuttaMethod(@getQ_nDxdt,cycle_count,dt,[y_a,y_a_tilde,zeros(3,1)],prm,isInterferrence);
        end

        prev_variable.y_a_tilde = y_a_tilde;

        %N
        y = N(cycle_count,dt,...
                [y_a_tilde,input(:,2),zeros(3,1)],...
                prm);
end