function nextStateVariable = cMoritaRungeKuttaMethod(obj,u)

    if ~obj.prm.settings.isRugekuttaMethodUse
        dxdt = obj.getDxdt(u,obj.stateVariable);
        obj.dxdtPreserved = dxdt; 
        nextStateVariable = obj.stateVariable + obj.dt*+dxdt;
    else
        k1 = obj.getDxdt(u,obj.stateVariable);
        obj.dxdtPreserved = k1;
        k2 = obj.getDxdt(u,obj.stateVariable + 0.5*obj.dt*k1);
        k3 = obj.getDxdt(u,obj.stateVariable + 0.5*obj.dt*k2);
        k4 = obj.getDxdt(u,obj.stateVariable + obj.dt*k3);
        nextStateVariable = obj.stateVariable + (k1 + 2*k2 + 2*k3 +k4) * (obj.dt/6);
    end

end