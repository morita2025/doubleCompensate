close all
 
% 
% %plotData
% tempData = transpose(variable.y([1,3],1:end));
% timeData = transpose(t);
% inputData = transpose(variable.u([1,3],:)); 
% RefData = transpose(variable.ref([1,3],:)); 
% plotTempData = [RefData, tempData];
% plotInputData = inputData;


%% plot


isExperimentGraph = true;
if isExperimentGraph % makeGraphの上書き
    load("C:\Users\mykot\OneDrive - Tokyo University of Agriculture and Technology\60MATLAB_sagyou\MicroreactorSystem2_morita\data\20250519stepOk.mat");
    tempData = data.temperature.sens(5:end,[4,6,1,3]);
    timeData = data.time(5:end,:);
    RefData = data.temperature.ref(4:end,:);
    inputData = data.current.controlInput(4:end,[1,3]);

    plotTempData = tempData;%[RefData, tempData];
    plotInputData = inputData;
    t=transpose(timeData);
else
    % plotTempData = prm.settings.outsideTemperature-[transpose(variable.ref([1,3],:)) ,transpose(variable.y([1,3],:))]; %
    % plotInputData = transpose(variable.u([1,3],:)); %
end


% makeGraph
FILE_IS_SAVE=true;
graphToolPath="C:\Users\mykot\OneDrive - Tokyo University of Agriculture and Technology\60MATLAB_sagyou\makeGraph";
addpath(graphToolPath);
DATA_DIR_PATH = "C:\Users\mykot\OneDrive - Tokyo University of Agriculture and Technology\40発表ゼミ\jisaku\figure\2025_5_19\"; %exp 
OUT_DIR_PATH = "C:\Users\mykot\OneDrive - Tokyo University of Agriculture and Technology\40発表ゼミ\jisaku\figure\2025_5_19\";
TEMPERATURE_GRAPH_TITLE = ["step_temperature_2"];
TEMPERATURE_LINE_NAME = ["$\mathrm{Part} \mathrm{A_1}$","$$\mathrm{Part} \mathrm{A_3}$","$\mathrm{Part} \mathrm{W_1}$","$$\mathrm{Part} \mathrm{W_3}$"];
TEMPERATURE_LINE_WIDTH = [2,2,2,2];
TEMPERATURE_LABEL_NAME = ["Time [$\mathrm{s}]$","Temperature [$^{\circ}\mathrm{C}]$"];
TEMPERATURE_LINE_STYLE = ["-","-","-","-",];
LINE_COLOR = transpose([[0.4 1 1]; [1 0.4 1]; [0 0 1]; [1 0 0]]);


CONTROLINPUT_GRAPH_TITLE = ["step_input"];
CONTROLINPUT_LINE_NAME = ["$u_1$","$u_3$"];
CONTROLINPUT_LABEL_NAME = ["Time [$\mathrm{s}]$","Current [$\mathrm{A}$]"];
%温度
makeGraph(t',plotTempData, ...
                    "lineName",TEMPERATURE_LINE_NAME, ...
                    "lineStyle",TEMPERATURE_LINE_STYLE, ...
                    "lineWidth",TEMPERATURE_LINE_WIDTH, ...
                    "labelName",TEMPERATURE_LABEL_NAME, ...
                    "graphName",TEMPERATURE_GRAPH_TITLE, ...
                    "location","north",...
                    "Color",LINE_COLOR,...
                    "lineWidth",[1,1,1,1], ..."yLimit",[19 27],...
                    "isSave",FILE_IS_SAVE,"outDir",OUT_DIR_PATH, ...
                    "fontSize",20,"LabelFontSize",30,"saveFileExt","png");
%制御入力
makeGraph(t',plotInputData , ...
                    "lineName",CONTROLINPUT_LINE_NAME, ...
                    "lineWidth",[1,1], ...
                    "labelName",CONTROLINPUT_LABEL_NAME, ...
                    "yLimit",[-0.2,1.8],...
                    "location","northwest",...
                    "graphName",CONTROLINPUT_GRAPH_TITLE, ...
                    "isSave",FILE_IS_SAVE,"outDir",OUT_DIR_PATH, ...
                    "fontSize",20,"LabelFontSize",30,"saveFileExt","png");


