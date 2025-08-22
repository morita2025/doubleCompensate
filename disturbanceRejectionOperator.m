function [f_1, f_2, prev_variable] = disturbanceRejectionOperator(cycleCount,dt,input,prev_variable,prm)
        y_a = input(:,1);
        u = input(:,2); %uは一ステップ前
        y = input(:,3);


        %D_inv
        y_a_tilde = D_inv(cycleCount,dt,[u,prev_variable.y_a_tilde_prev,zeros(3,1)],prm,false);
        d1_hat = y_a - y_a_tilde ;%d1の推定

        %g(NQ(d2)の推定)
        y_w_tilde = N(cycleCount,dt,[prev_variable.y_a,prev_variable.y_w_tilde_prev,zeros(3,1)],prm,false);
        g = y - y_w_tilde;

        %M
        plantStateVariableM = prev_variable.invPlantStateVariableM;
           if prm.settings.isRugekuttaMethodUse == 0
                dxdt2 = -prm.lowPassFilterTimePrm.p * plantStateVariableM(:,1) + prm.lowPassFilterTimePrm.p * d1_hat;
           else
                k1 = -prm.lowPassFilterTimePrm.p * plantStateVariableM(:,1) + prm.lowPassFilterTimePrm.p * d1_hat;
                k2 = -prm.lowPassFilterTimePrm.p * (0.5*dt*k1 + plantStateVariableM(:,1)) + prm.lowPassFilterTimePrm.p * d1_hat;
                k3 = -prm.lowPassFilterTimePrm.p * (0.5*dt*k2 + plantStateVariableM(:,1)) + prm.lowPassFilterTimePrm.p * d1_hat;
                k4 = -prm.lowPassFilterTimePrm.p * (dt*k3 + plantStateVariableM(:,1)) + prm.lowPassFilterTimePrm.p * d1_hat;
                dxdt2 = (k1 + 2*k2 + 2*k3 +k4) /6;
           end


        f_m = prm.lowPassFilterTimePrm.D.* dt*dxdt2 +plantStateVariableM(:,1);
        plantStateVariableM(:,1) = dxdt2*dt + plantStateVariableM(:,1);
        f_m = f_m .* prm.MOperatorConstPrm;
        plantStateVariableM(:,2) = f_m;

        %N
        f1Noutput = N(cycleCount,dt,[d1_hat,prev_variable.y_w_prev,zeros(3,1)],prm,0);
        if prm.settings.isD2Compensate
            f_a = f1Noutput + g;
        else
            f_a = f1Noutput;
        end

        % plantStateVariable = prev_variable.invPlantStateVariable;
        % p_3 = prm.A_w;
        % p_4 = prm.interferrenceConstPrm.y_aw ./ [prm.mw_cw13; prm.mw_cw2; prm.mw_cw13];
        % dxdt = [-2*prm.lowPassFilterTimePrm.p*plantStateVariable(:,1) - (prm.lowPassFilterTimePrm.p^2)*plantStateVariable(:,2) + (prm.lowPassFilterTimePrm.p^2)*g, plantStateVariable(:,1)];
        % plantStateVariable = dxdt*dt + plantStateVariable;
        % 
        % 
        % a1 = (prm.lowPassFilterTimePrm.N ./ p_4);
        % a2 = ((1./prm.lowPassFilterTimePrm.N) + p_3)  .*plantStateVariable(:,1);
        % a3 = (p_3./prm.lowPassFilterTimePrm.N).*plantStateVariable(:,2);
        % d_2hat =   a1.* ( dxdt(:,1) + (a2 + a3)); %d_2hatも合ってる


        %修正
        % plantStateVariable = prev_variable.invPlantStateVariable;
        % p_3 = prm.A_w;
        % p_4 = prm.interferrenceConstPrm.y_aw ./ [prm.mw_cw13; prm.mw_cw2; prm.mw_cw13];
        % dxdt = -prm.lowPassFilterTimePrm.p * plantStateVariable(:,1) + prm.lowPassFilterTimePrm.p * g;
        % plantStateVariable(:,1) = dxdt*dt + plantStateVariable(:,1);
        % d_2hat = 1./p_4 .* (p_3.*plantStateVariable(:,1) + dxdt);
        % 


        %修正0602 (NQ)^{-1} F 0707ルンゲクッタ法に対応
        plantStateVariable = prev_variable.invPlantStateVariable;
        p_3 = prm.A_w;
        p_4 = prm.interferrenceConstPrm.y_aw ./ [prm.mw_cw13; prm.mw_cw2; prm.mw_cw13];
        %F2
        % dx2dt = -2*prm.lowPassFilterTimePrm.p2 * plantStateVariable(:,2) - prm.lowPassFilterTimePrm.p2^2 * plantStateVariable(:,1) + prm.lowPassFilterTimePrm.p2^2 * g;
           if prm.settings.isRugekuttaMethodUse == 0
                dx2dt = -2*prm.lowPassFilterTimePrm.p2 * plantStateVariable(:,2) - prm.lowPassFilterTimePrm.p2^2 * plantStateVariable(:,1) + prm.lowPassFilterTimePrm.p2^2 * g;
                %integral
                plantStateVariable(:,2) = dx2dt*dt + plantStateVariable(:,2);
                plantStateVariable(:,1) = plantStateVariable(:,2)*dt + plantStateVariable(:,1);
           else
                k1_1 = -2*prm.lowPassFilterTimePrm.p2 * plantStateVariable(:,2) - prm.lowPassFilterTimePrm.p2^2 * plantStateVariable(:,1) + prm.lowPassFilterTimePrm.p2^2 * g;
                k1_2 = plantStateVariable(:,2);
                k2_1 = -2*prm.lowPassFilterTimePrm.p2 * (0.5*dt*k1_1 + plantStateVariable(:,2)) - prm.lowPassFilterTimePrm.p2^2 * (0.5*dt*k1_2 +plantStateVariable(:,1)) + prm.lowPassFilterTimePrm.p2^2 * g;
                k2_2 = 0.5*dt*k1_1+plantStateVariable(:,2);
                k3_1 = -2*prm.lowPassFilterTimePrm.p2 * (0.5*dt*k2_1 + plantStateVariable(:,2)) - prm.lowPassFilterTimePrm.p2^2 * (0.5*dt*k2_2 +plantStateVariable(:,1)) + prm.lowPassFilterTimePrm.p2^2 * g;
                k3_2 = 0.5*dt*k2_1+plantStateVariable(:,2);
                k4_1 = -2*prm.lowPassFilterTimePrm.p2 * (dt*k3_1 + plantStateVariable(:,2)) - prm.lowPassFilterTimePrm.p2^2 * (dt*k3_2 +plantStateVariable(:,1)) + prm.lowPassFilterTimePrm.p2^2 * g;
                k4_2 = dt*k3_1+plantStateVariable(:,2);
                dx2dt_1 = (k1_1 + 2*k2_1 + 2*k3_1 +k4_1) /6;
                dx2dt_2 = (k1_2 + 2*k2_2 + 2*k3_2 +k4_2) /6;
                %integral
                plantStateVariable(:,2) = dx2dt_1*dt + plantStateVariable(:,2);
                plantStateVariable(:,1) = dx2dt_2*dt + plantStateVariable(:,1);
                dx2dt = dx2dt_1;
           end

        %(NQ)^{-1}
        d_2hat = 1./p_4 .* ( prm.lowPassFilterTimePrm.D.*dx2dt + (prm.lowPassFilterTimePrm.D.*p_3 + 1).*plantStateVariable(:,2) + p_3.* plantStateVariable(:,1));




        prev_variable.y_w_prev = f1Noutput; %f1のN
        prev_variable.y_a = y_a;
        prev_variable.g = g;
        prev_variable.y_w_tilde_prev = y_w_tilde;
        prev_variable.y_a_tilde_prev = y_a_tilde;
        prev_variable.invPlantStateVariable = plantStateVariable;
        prev_variable.invPlantStateVariableM = plantStateVariableM;
        prev_variable.d_prev = d1_hat;%使う
        prev_variable.debug = d_2hat;
        
        if prm.settings.isD1Compensate
            f_1 = f_a -f_m;
        else
            f_1=0;
        end

        if prm.settings.isD2Compensate
            f_2 = d_2hat.*prm.MOperatorConstPrm;
        else
            f_2=0;
        end


        

end