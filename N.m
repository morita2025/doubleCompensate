function y = N(cycle_count,dt,prev_variable,prm,isInterferrence)
         arguments
             cycle_count
             dt
             prev_variable
             prm
             isInterferrence=prm.settings.isInterferrence
         end


        if prm.settings.isRugekuttaMethodUse == 0
            y = moritaEulerMethod(@getNDxdt,cycle_count,dt,prev_variable,prm,isInterferrence);
        end

        if prm.settings.isRugekuttaMethodUse == 1
            y = moritaRungekuttaMethod(@getNDxdt,cycle_count,dt,prev_variable,prm,isInterferrence);
        end

end






