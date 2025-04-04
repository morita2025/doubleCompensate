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
        f_m = d1_hat .* prm.MOperatorConstPrm;

        %N
        f1Noutput = N(cycleCount,dt,[d1_hat,prev_variable.y_w_prev,zeros(3,1)],prm,0);
        if prm.settings.isD2Compensate
            f_a = f1Noutput + g;
        else
            f_a = f1Noutput;
        end

        plantStateVariable = prev_variable.invPlantStateVariable;
        p_3 = prm.A_w;
        p_4 = prm.interferrenceConstPrm.y_aw ./ [prm.mw_cw13; prm.mw_cw2; prm.mw_cw13];
        dxdt = [-2*prm.lowPassFilterTimePrm.p*plantStateVariable(:,1) - (prm.lowPassFilterTimePrm.p^2)*plantStateVariable(:,2) + (prm.lowPassFilterTimePrm.p^2)*g, plantStateVariable(:,1)];
        plantStateVariable = dxdt*dt + plantStateVariable;


        a1 = (prm.lowPassFilterTimePrm.N ./ p_4);
        a2 = ((1./prm.lowPassFilterTimePrm.N) + p_3)  .*plantStateVariable(:,1);
        a3 = (p_3./prm.lowPassFilterTimePrm.N).*plantStateVariable(:,2);
        d_2hat =   a1.* ( dxdt(:,1) + (a2 + a3)); %d_2hatも合ってる

        prev_variable.y_w_prev = f1Noutput; %f1のN
        prev_variable.y_a = y_a;
        prev_variable.y_w_tilde_prev = y_w_tilde;
        prev_variable.y_a_tilde_prev = y_a_tilde;
        prev_variable.invPlantStateVariable = plantStateVariable;
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