function dxdt = getQ_nDxdt(cycleCount,dt,input,prm,isInterferrence)
         arguments
             cycleCount
             dt
             input
             prm
             isInterferrence=prm.settings.isInterferrence
         end


        x_1 = input(:,1);
        x_2 = input(:,2);

        dxdt = (x_1 - x_2) / prm.lowPassFilterTimePrm.N;



end