function next_cycle_variable = moritaEulerMethod(functionPointer,cycle_count,dt,prev_variable,prm,isInterferrence)
         arguments
             functionPointer
             cycle_count
             dt
             prev_variable
             prm
             isInterferrence=prm.settings.isInterferrence
         end

    dxdt = functionPointer(cycle_count,dt,prev_variable,prm);
    next_cycle_variable = dxdt * dt + prev_variable(:,2);


end