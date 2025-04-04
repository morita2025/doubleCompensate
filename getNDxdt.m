function dxdt = getNDxdt(cycle_count,dt,input,prm,isInterferrence)
         arguments
             cycle_count
             dt
             input
             prm
             isInterferrence=prm.settings.isInterferrence
         end
            y_a = input(:,1);
            y_prev = input(:,2);
            %入力熱量
            dxdt = prm.interferrenceConstPrm.y_aw.*y_a./[prm.mw_cw13; prm.mw_cw2; prm.mw_cw13];
            %プラント
            dxdt = dxdt - prm.A_w.*y_prev;
            %干渉
            if isInterferrence
                dxdt = dxdt + prm.interferrenceConstPrm.y_w.*[y_prev(2); y_prev(1)+y_prev(3); y_prev(2)]./[prm.mw_cw13; prm.mw_cw2; prm.mw_cw13];
            end
        
end