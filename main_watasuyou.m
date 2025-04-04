clear ;
close all;
max_time = 100;
dt = 1;
t = 0:dt:max_time;
ref = [3.8; 4];

prm =  CalcOperatorPrm_kato(outsideTemperature=28,max_time=max_time,i_max=2,i_min=0,heatTransferCoef_water=270,...
                           isRugekuttaMethodUse=1,isInterferrence=false);
variable= getVariableFunction(length(t),ref);
operatorTempVariable = struct("B_inv",struct("y_w_prev",zeros(3,1),"x_1_prev",zeros(3,1),"x_2_prev",zeros(3,1),"x_3_prev",zeros(3,1),"x_debug_prev",zeros(3,1),"x_debug2_prev",zeros(3,1),"x_debug3_prev",zeros(3,1)),...
                              "N_tilde",struct("y_a_tilde",zeros(3,1)),...
                              "D_tilde_inv",struct( ),...
                              "disturbanceRejectionOperator",struct("y_w_prev",zeros(3,1),"y_a_tilde_prev",zeros(3,1),"d_c_prev",zeros(3,1),"d_prev",zeros(3,1)));


for cycleCount = 1:length(t)
    if cycleCount > 1 
        %% プラント
        variable.y_a(:,cycleCount)= D_inv(cycleCount,dt,...
                        [variable.u(:,cycleCount-1),variable.y_a(:,cycleCount-1),variable.y(:,cycleCount-1)],...
                        prm);

        variable.y(:,cycleCount) = N(cycleCount,dt,...
                [variable.y_a(:,cycleCount-1),variable.y(:,cycleCount-1),zeros(3,1)],...
                prm);



    end
end


