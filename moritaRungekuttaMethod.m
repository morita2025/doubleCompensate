function next_cycle_variable = moritaRungekuttaMethod(functionPointer,cycle_count,dt,prev_variable,prm,isInterferrence)
         arguments
             functionPointer
             cycle_count
             dt
             prev_variable
             prm
             isInterferrence=prm.settings.isInterferrence
         end

        half_dt = dt/2;
        k1 = functionPointer(cycle_count,dt,prev_variable,prm,isInterferrence);
        k2 = functionPointer(cycle_count,dt,prev_variable+half_dt*[zeros(3,1),k1,zeros(3,1)],prm,isInterferrence);
        k3 = functionPointer(cycle_count,dt,prev_variable+half_dt*[zeros(3,1),k2,zeros(3,1)],prm,isInterferrence);
        k4 = functionPointer(cycle_count,dt,prev_variable+dt*[zeros(3,1),k3,zeros(3,1)],prm,isInterferrence);
    
        next_cycle_variable = (k1 + 2*k2 + 2*k3 +k4) * (dt/6) + prev_variable(:,2);


end