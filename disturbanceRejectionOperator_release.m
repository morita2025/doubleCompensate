function [f_1, f_2, prev_variable] = disturbanceRejectionOperator_release(cycleCount,dt,input,prev_variable,prm)
        y_a = input(:,1);
        u = input(:,2); %uは一ステップ前
        y = input(:,3);


        %D_inv
        y_a_tilde = D_inv(cycleCount,dt,[u,prev_variable.y_a_tilde_prev,zeros(3,1)],prm,false);
        d1_hat = y_a - y_a_tilde ;%d1の推定

        %g(NQ(d2)の推定)
        y_w_tilde = N(cycleCount,dt,[prev_variable.y_a,prev_variable.y_w_tilde_prev,zeros(3,1)],prm,false);
        g = y - y_w_tilde;

        %20250826　新M (F1が2次遅れ)
        plantStateVariableInvQF1 = zeros(3,1);
        plantStateVariableInvQF1([1,3]) = prev_variable.instance.controllerinvQF1.calcNextCycle(d1_hat([1,3]));
        f_m = prm.MOperatorConstPrm .* plantStateVariableInvQF1;

        %N
        f1Noutput = N(cycleCount,dt,[d1_hat,prev_variable.y_w_prev,zeros(3,1)],prm,0);


        %20250821追加
        Aoutput=zeros(3,1);
        if prm.settings.isD2Compensate
            Aoutput([1,3]) = prev_variable.instance.controllerA.calcNextCycle(f1Noutput([1,3])+ g([1,3]));
        else
            Aoutput([1,3]) = prev_variable.instance.controllerA.calcNextCycle(f1Noutput([1,3]));
        end
        f_a = Aoutput;

        %20250826F2を三次遅れに変更
        d_2hat=zeros(3,1);
        d_2hat([1,3]) = prev_variable.instance.controllerinvTildeNF2.calcNextCycle(g([1,3]));



        prev_variable.y_w_prev = f1Noutput; %f1のN
        prev_variable.y_a = y_a;
        prev_variable.g = g;
        prev_variable.y_w_tilde_prev = y_w_tilde;
        prev_variable.y_a_tilde_prev = y_a_tilde;
        prev_variable.d_prev = d1_hat;%使う
        prev_variable.debug = d_2hat;
        prev_variable.Aoutput_prev = Aoutput;
        
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