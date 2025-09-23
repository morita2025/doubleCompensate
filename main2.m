clear ;
close all;
max_time = 1200;
dt = 1;
t = 0:dt:max_time;
ref = [4; 4.5];
addpath(pwd + "\Abstract\")

prm =  CalcOperatorPrm_kato(outsideTemperature=28,max_time=max_time,i_max=2,i_min=0,heatTransferCoef_water=270,tau=30,p=0.2,p2=0.1,p_A=1,...
                           isRugekuttaMethodUse=1,isInterferrence=true,isD1Compensate=true,isD2Compensate=true);
variable= getVariableFunction(length(t),ref);
operatorTempVariable = struct("B_inv",struct("y_w_prev",zeros(3,1),"x_1_prev",zeros(3,1),"x_2_prev",zeros(3,1),"Aoutput_prev",zeros(3,1),...
                                            "x_3_prev",zeros(3,1),"x_debug_prev",zeros(3,1),"x_debug2_prev",zeros(3,1),"x_debug3_prev",zeros(3,1),...
                                            "instance",struct("controllerA",ControllerA(prm=prm,dt=dt,cycleNum=length(t)))),...
                          "N_tilde",struct("y_a_tilde",zeros(3,1)),...
                          "D_tilde_inv",struct( ),...
                          "disturbanceRejectionOperator",struct("y_w_prev",zeros(3,1),"y_a_tilde_prev",zeros(3,1),"y_a",zeros(3,1),...
                                                                "invPlantStateVariable",zeros(3,2),"invPlantStateVariableM",zeros(3,2),"debug",zeros(3,2),...
                                                                "g",zeros(3,1),"Aoutput_prev",zeros(3,1),...
                                                                "instance",struct("controllerA",ControllerA(prm=prm,dt=dt,cycleNum=length(t)),"controllerinvQF1",ControllerinvQF1(prm=prm,dt=dt,cycleNum=length(t)),...
                                                                "controllerinvTildeNF2",ControllerinvTildeNF2(prm=prm,dt=dt,cycleNum=length(t))),...
                                                                "y_w_tilde_prev",zeros(3,1),"d_c_prev",zeros(3,1),"d_prev",zeros(3,1)));

refTimePrm = 1/10;
refChangeTime = 400;
variable.ref(:,:) = [ref(1); 0; ref(2)].* (1 - exp(-refTimePrm*t)  );
% variable.ref(:,refChangeTime:end)  = [5; 0; 5.5] .* ones(3,max_time-refChangeTime+2);

% variable.tubeGairan([1,3],400:end) = -0.5;
variable.tubeGairan([1,3],800:end) = -1;
variable.almiGairan([1,3],450:end) = -0.5;
% variable.almiGairan([1,3],800:end) = -1;



%Aインスタンス実験
instanceInvQF1 = ControllerinvQF1(prm=prm,dt=dt,cycleNum=length(t)); 
instanceInvTildeNF2 = ControllerinvTildeNF2(prm=prm,dt=dt,cycleNum=length(t)); 
instanceA = ControllerA(prm=prm,dt=dt,cycleNum=length(t)); 





for cycleCount = 1:length(t)
    if cycleCount > 0 

        % 外乱除去制御系
        if cycleCount >1
            [variable.f_1(:,cycleCount), variable.f_2(:,cycleCount), operatorTempVariable.disturbanceRejectionOperator] = disturbanceRejectionOperator_release(cycleCount,dt,...
                    [variable.y_a(:,cycleCount),variable.u(:,cycleCount-1),variable.y(:,cycleCount)],operatorTempVariable.disturbanceRejectionOperator,...
                    prm);
        end

        %% 制御器
        variable.r_01(:,cycleCount) = variable.ref(:,cycleCount) - variable.f_2(:,cycleCount);
        variable.r_02(:,cycleCount) = variable.r_01(:,cycleCount) + variable.f_1(:,cycleCount);

        variable.b([1,3],cycleCount) =  instanceA.calcNextCycle(variable.y([1,3],cycleCount));
        variable.e(:,cycleCount) =  variable.r_02(:,cycleCount) -variable.b(:,cycleCount);  


        [variable.u(:,cycleCount), operatorTempVariable.B_inv] = B_inv(cycleCount,dt,...
                variable.e(:,cycleCount),operatorTempVariable.B_inv,...
                prm);



        %% プラント
        if cycleCount ~= length(t)
            variable.y_a_asterisk(:,cycleCount+1)= D_inv(cycleCount,dt,...
                            [variable.u(:,cycleCount),variable.y_a_asterisk(:,cycleCount),variable.y(:,cycleCount)],...
                            prm);
            %外乱(アルミ)
            variable.y_a(:,cycleCount+1) = variable.y_a_asterisk(:,cycleCount+1) + variable.almiGairan(:,cycleCount);

            %
            variable.y_asterisk(:,cycleCount+1) = N(cycleCount,dt,...
                    [variable.y_a(:,cycleCount),variable.y_asterisk(:,cycleCount),zeros(3,1)],...
                    prm);

            %外乱(チューブ)
            variable.y(:,cycleCount+1) = variable.y_asterisk(:,cycleCount+1) + variable.tubeGairan(:,cycleCount+1);
        
        end

        %% プラント(右分解)
        % [variable.y_f(:,cycleCount), operatorTempVariable.D_tilde_inv] = D_tilde_inv(cycleCount,dt,...
        %                 [variable.u(:,cycleCount-1),variable.y_f(:,cycleCount-1),variable.y(:,cycleCount-1)],operatorTempVariable.D_tilde_inv,...
        %                 prm,0);
        % 
        % [variable.y(:,cycleCount), operatorTempVariable.N_tilde] = N_tilde(cycleCount,dt,...
        %         [variable.y_f(:,cycleCount-1),variable.y(:,cycleCount-1),zeros(3,1)],operatorTempVariable.N_tilde,...
        %         prm,0);

        %% NM^{-1}
        % u_ast = M(cycleCount,dt,variable.ref(:,cycleCount),prm);
        % variable.y_f(:,cycleCount+1) = N(cycleCount,dt,...
        %         [u_ast,variable.y_f(:,cycleCount),zeros(3,1)],...
        %         prm);

        %% debug
        % y_a_dot = getD_invDxdt(cycleCount,dt,...
        %                 [variable.u(:,cycleCount),variable.y_a(:,cycleCount),variable.y(:,cycleCount)],...
        %                 prm,0);
        variable.y_f(:,cycleCount) = operatorTempVariable.disturbanceRejectionOperator.d_prev(:,1);
        variable.y_g(:,cycleCount) = operatorTempVariable.disturbanceRejectionOperator.invPlantStateVariable(:,1);
        variable.y_v(:,cycleCount) = operatorTempVariable.disturbanceRejectionOperator.debug(:,1);
        variable.y_h(:,cycleCount) = operatorTempVariable.disturbanceRejectionOperator.g(:,1);

    end



    %Aの実験
    % Aoutput(:,cycleCount) = instanceInvQF1.calcNextCycle([1;1]);
    Aoutput(:,cycleCount) = instanceInvTildeNF2.calcNextCycle([1;1]);
end

% plot(Aoutput(2,:))



% plot(variable.y_g(1,:),"Color","b","LineWidth",2)
% hold on
% plot(variable.y_v(1,:))


% figure()
% plot(t,variable.r_01(1,:),"r")%almi
% hold on
% plot(t,variable.ref(1,1:end),"b") %tube

% figure(2)
% plot(t,variable.f_1(1,:) - variable.f_2(1,:))
%ブレークポイントに関するショートカットを調べる


%リポジトリ, init, add, commit(＋初回呪文), push, bramch, checkout, pull





%plotData
tempData = transpose(variable.y([1,3],1:end));
timeData = transpose(t);
inputData = transpose(variable.u([1,3],:)); 
RefData = transpose(variable.ref([1,3],:)); 
plotTempData = [RefData, tempData];
plotInputData = inputData;



%% plot


isExperimentGraph = false;
if isExperimentGraph % makeGraphの上書き
    load("C:\Users\mykot\OneDrive - Tokyo University of Agriculture and Technology (1)\60MATLAB_sagyou\MicroreactorSystem2_morita\data\20250714_kato_5.mat");
    tempData = data.temperature.sens(5:end,[1,3]);
    timeData = data.time(5:end,:);
    RefData = data.temperature.ref(4:end,:);
    inputData = data.current.controlInput(4:end,[1,3]);

    plotTempData = [RefData, tempData];
    plotInputData = inputData;
    t=transpose(timeData);
else
    plotTempData = prm.settings.outsideTemperature-[transpose(variable.ref([1,3],:)) ,transpose(variable.y([1,3],:))]; %
    plotInputData = transpose(variable.u([1,3],:)); %
end


% makeGraph
FILE_IS_SAVE=false;
graphToolPath="C:\Users\mykot\OneDrive - Tokyo University of Agriculture and Technology\60MATLAB_sagyou\makeGraph";
addpath(graphToolPath);
DATA_DIR_PATH = "C:\Users\mykot\OneDrive - Tokyo University of Agriculture and Technology\40発表ゼミ\jisaku\figure\2025_5_19\"; %exp 
OUT_DIR_PATH = "C:\Users\mykot\OneDrive - Tokyo University of Agriculture and Technology\40発表ゼミ\jisaku\figure\2025_7_20\";
TEMPERATURE_GRAPH_TITLE = ["gairan"];
TEMPERATURE_LINE_NAME = ["$T_{0}-r_1$","$T_{0}-r_3$","$\mathrm{Part} \mathrm{W_1}$","$$\mathrm{Part} \mathrm{W_3}$"];
TEMPERATURE_LINE_WIDTH = [2,2,2,2];
TEMPERATURE_LABEL_NAME = ["Time [$\mathrm{s}]$","Temperature [$^{\circ}\mathrm{C}]$"];
TEMPERATURE_LINE_STYLE = ["--","--","-","-",];


CONTROLINPUT_GRAPH_TITLE = ["temperature_kato"];
CONTROLINPUT_LINE_NAME = ["$u_1$","$u_3$"];
CONTROLINPUT_LABEL_NAME = ["Time [$\mathrm{s}]$","Current [$\mathrm{A}$]"];
%温度
makeGraph(t',plotTempData, ...
                    "lineName",TEMPERATURE_LINE_NAME, ...
                    "lineStyle",TEMPERATURE_LINE_STYLE, ...
                    "lineWidth",TEMPERATURE_LINE_WIDTH, ...
                    "labelName",TEMPERATURE_LABEL_NAME, ...
                    "graphName",TEMPERATURE_GRAPH_TITLE, ...
                    "location","northeast",...
                    "lineWidth",[2,2,2,2], ...
                    "yLimit",[22 28],...
                    "isSave",FILE_IS_SAVE,"outDir",OUT_DIR_PATH, ...
                    "fontSize",20,"LabelFontSize",30,"saveFileExt","png");
% %制御入力
% makeGraph(t',plotInputData , ...
%                     "lineName",CONTROLINPUT_LINE_NAME, ...
%                     "lineWidth",[1,1], ...
%                     "labelName",CONTROLINPUT_LABEL_NAME, ...
%                     "yLimit",[-0.2,1.4],...
%                     "location","northwest",...
%                     "graphName",CONTROLINPUT_GRAPH_TITLE, ...
%                     "isSave",FILE_IS_SAVE,"outDir",OUT_DIR_PATH, ...
%                     "fontSize",20,"LabelFontSize",30,"saveFileExt","png");


