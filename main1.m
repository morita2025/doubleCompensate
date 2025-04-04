clear ;
close all;
max_time = 1000;
dt = 1;
t = 0:dt:max_time;
ref = [3.8; 4];

prm =  CalcOperatorPrm_kato(outsideTemperature=28,max_time=max_time,i_max=2,i_min=0,heatTransferCoef_water=270,p=0.6,...
                           isRugekuttaMethodUse=1,isInterferrence=true,isD1Compensate=true,isD2Compensate=true);
variable= getVariableFunction(length(t),ref);
operatorTempVariable = struct("B_inv",struct("y_w_prev",zeros(3,1),"x_1_prev",zeros(3,1),"x_2_prev",zeros(3,1),"x_3_prev",zeros(3,1),"x_debug_prev",zeros(3,1),"x_debug2_prev",zeros(3,1),"x_debug3_prev",zeros(3,1)),...
                          "N_tilde",struct("y_a_tilde",zeros(3,1)),...
                          "D_tilde_inv",struct( ),...
                          "disturbanceRejectionOperator",struct("y_w_prev",zeros(3,1),"y_a_tilde_prev",zeros(3,1),"y_a",zeros(3,1),...
                                                                "invPlantStateVariable",zeros(3,2),"debug",zeros(3,2),...
                                                                "y_w_tilde_prev",zeros(3,1),"d_c_prev",zeros(3,1),"d_prev",zeros(3,1)));

refTimePrm = 1/10;
refChangeTime = 400;
variable.ref(:,:) = [ref(1); 0; ref(2)].* (1 - exp(-refTimePrm*t)  );
variable.ref(:,refChangeTime:end)  = [5; 0; 5.5] .* ones(3,max_time-refChangeTime+2);

variable.tubeGairan([1,3],600:end) = -1;%*repmat(sin(0.01*t(1:402)),2,1);
variable.almiGairan([1,3],800:end) = -1*ones(2,202);


for cycleCount = 1:length(t)
    if cycleCount > 0 

        % 外乱除去制御系
        if cycleCount >1
            [variable.f_1(:,cycleCount), variable.f_2(:,cycleCount), operatorTempVariable.disturbanceRejectionOperator] = disturbanceRejectionOperator(cycleCount,dt,...
                    [variable.y_a(:,cycleCount),variable.u(:,cycleCount-1),variable.y(:,cycleCount)],operatorTempVariable.disturbanceRejectionOperator,...
                    prm);
        end

        %% 制御器
        variable.r_01(:,cycleCount) = variable.ref(:,cycleCount) - variable.f_2(:,cycleCount);
        variable.r_02(:,cycleCount) = variable.r_01(:,cycleCount) + variable.f_1(:,cycleCount);

        variable.b(:,cycleCount) =  variable.y(:,cycleCount);
        variable.e(:,cycleCount) =  variable.r_02(:,cycleCount) -variable.b(:,cycleCount);  


        [variable.u(:,cycleCount), operatorTempVariable.B_inv] = B_inv(cycleCount,dt,...
                variable.e(:,cycleCount),operatorTempVariable.B_inv,...
                prm);



        %% プラント
        variable.y_a_asterisk(:,cycleCount+1)= D_inv(cycleCount,dt,...
                        [variable.u(:,cycleCount),variable.y_a_asterisk(:,cycleCount),variable.y(:,cycleCount)],...
                        prm);
        %外乱(アルミ)
        variable.y_a(:,cycleCount+1) = variable.y_a_asterisk(:,cycleCount+1) + variable.almiGairan(:,cycleCount);
        %外乱(チューブ)
        variable.ya_gairan(:,cycleCount) = variable.y_a(:,cycleCount) + variable.tubeGairan(:,cycleCount);

        %
        variable.y(:,cycleCount+1) = N(cycleCount,dt,...
                [variable.ya_gairan(:,cycleCount),variable.y(:,cycleCount),zeros(3,1)],...
                prm);

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

    end
end
figure()
plot(t,variable.r_01(1,:),"r")%almi
hold on
plot(t,variable.ref(1,1:end),"b") %tube

% figure(2)
% plot(t,variable.f_1(1,:) - variable.f_2(1,:))
%ブレークポイントに関するショートカットを調べる


%リポジトリ, init, add, commit(＋初回呪文), push, bramch, checkout, pull





%plotData
tempData = transpose(variable.y([1,3],1:end-1));
timeData = transpose(t);
inputData = transpose(variable.u([1,3],:)); 
RefData = transpose(variable.ref([1,3],:)); 
plotTempData = [RefData, tempData];
plotInputData = inputData;



% 
% %% makeGraphの上書き
% load("C:\Users\mykot\OneDrive - Tokyo University of Agriculture and Technology\60MATLAB_sagyou\MicroreactorSystem1\data\0707_hikaku.mat");
% tempData = data.data.temperature.sensor(3:end,[1,3]);
% timeData = data.data.time(3:end,:);
% RefData = data.data.temperature.ref(3:end,:);
% inputData = data.data.other(3:end,[1,2]);
% 
% plotTempData = [RefData, tempData];
% plotInputData = inputData;
% t=transpose(timeData);
% 
% 
% 


% makeGraph
FILE_IS_SAVE=false;
graphToolPath="C:\Users\mykot\OneDrive - Tokyo University of Agriculture and Technology\60MATLAB_sagyou\makeGraph";
addpath(graphToolPath);
DATA_DIR_PATH = "C:\Users\mykot\OneDrive - Tokyo University of Agriculture and Technology\60MATLAB_sagyou\kato_simulation\graph"; %exp 
OUT_DIR_PATH = "C:\Users\mykot\OneDrive - Tokyo University of Agriculture and Technology\60MATLAB_sagyou\kato_simulation\graph\";
TEMPERATURE_GRAPH_TITLE = ["temperature(ronbun_hikaku)"];
TEMPERATURE_LINE_NAME = ["$T_{0_2}-r^{\ast}_1$","$T_{0_2}-r^{\ast}_3$","$\mathrm{Part} \mathrm{W_1}$","$$\mathrm{Part} \mathrm{W_3}$"];
TEMPERATURE_LINE_WIDTH = [1,1,2,2];
TEMPERATURE_LABEL_NAME = ["Time [$\mathrm{s}]$","Temperature [$^{\circ}\mathrm{C}]$"];
TEMPERATURE_LINE_STYLE = ["--","--","-","-",];


CONTROLINPUT_GRAPH_TITLE = ["controlInput(ronbun_hikaku)"];
CONTROLINPUT_LINE_NAME = ["$u_1$","$u_3$"];
CONTROLINPUT_LABEL_NAME = ["Time [$\mathrm{s}]$","Current [$\mathrm{A}$]"];
%温度
% plotTempData = prm.settings.outsideTemperature-[transpose(variable.ref([1,3],:)) ,transpose(variable.y([1,3],:))]; %
makeGraph(t',plotTempData, ...
                    "lineName",TEMPERATURE_LINE_NAME, ...
                    "lineStyle",TEMPERATURE_LINE_STYLE, ...
                    "lineWidth",TEMPERATURE_LINE_WIDTH, ...
                    "labelName",TEMPERATURE_LABEL_NAME, ...
                    "graphName",TEMPERATURE_GRAPH_TITLE, ...
                    "location","southwest",...
                    "lineWidth",[1,1,1,1], ..."yLimit",[22.4 27.5],...
                    "isSave",FILE_IS_SAVE,"outDir",OUT_DIR_PATH, ...
                    "fontSize",20,"LabelFontSize",30,"saveFileExt","jpg");
%制御入力
% plotInputData = transpose(variable.u([1,3],:)); %
makeGraph(t',plotInputData , ...
                    "lineName",CONTROLINPUT_LINE_NAME, ...
                    "lineWidth",[1,1], ...
                    "labelName",CONTROLINPUT_LABEL_NAME, ..."yLimit",[-0.8,1.6],...
                    "location","northwest",...
                    "graphName",CONTROLINPUT_GRAPH_TITLE, ...
                    "isSave",FILE_IS_SAVE,"outDir",OUT_DIR_PATH, ...
                    "fontSize",20,"LabelFontSize",30,"saveFileExt","jpg");


