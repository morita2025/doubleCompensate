%% makeGraphの上書き 実機実験の結果をプロットする関数
% load("C:\Users\mykot\OneDrive - Tokyo University of Agriculture and Technology (1)\60MATLAB_sagyou\MicroreactorSystem1\data\0707_hikaku.mat");
load("C:\Users\mykot\OneDrive - Tokyo University of Agriculture and Technology (1)\60MATLAB_sagyou\kato_simulation\0703_1.mat");
tempData = data.data.temperature.sensor(3:end,[1,3]);
timeData = data.data.time(3:end,:);
RefData = data.data.temperature.ref(3:end,:);
inputData = data.data.current.ref(3:end,[1,3]);

timeDelayData = data.data.timeDelay(3:end,:);
execTimeData = data.data.UserProgramExecTime(3:end,:);

plotTempData = [RefData, tempData];
plotInputData = inputData;
t=transpose(timeData);





%% makeGraph
FILE_IS_SAVE=true;
graphToolPath="C:\Users\mykot\OneDrive - Tokyo University of Agriculture and Technology (1)\60MATLAB_sagyou\makeGraph";
addpath(graphToolPath);
DATA_DIR_PATH = "C:\Users\mykot\OneDrive - Tokyo University of Agriculture and Technology (1)\60MATLAB_sagyou\kato_simulation\graph"; %exp 
OUT_DIR_PATH = "C:\Users\mykot\OneDrive - Tokyo University of Agriculture and Technology (1)\60MATLAB_sagyou\kato_simulation\graph\";
% TEMPERATURE_GRAPH_TITLE = ["temperature"];
% TEMPERATURE_LINE_NAME = ["$T_{0}-r_1$","$T_{0}-r_3$","$\mathrm{Part} \mathrm{W_1}$","$$\mathrm{Part} \mathrm{W_3}$"];
% TEMPERATURE_LINE_WIDTH = [1,1,2,2];
% TEMPERATURE_LABEL_NAME = ["Time [$\mathrm{s}]$","Temperature [$^{\circ}\mathrm{C}]$"];
% TEMPERATURE_LINE_STYLE = ["--","--","-","-",];
% 
% 
% CONTROLINPUT_GRAPH_TITLE = ["controlInput"];
% CONTROLINPUT_LINE_NAME = ["$i_1$","$i_3$"];
% CONTROLINPUT_LABEL_NAME = ["Time [$\mathrm{s}]$","Current [$\mathrm{A}$]"];
% %温度
% makeGraph(t',plotTempData, ...
%                     "lineName",TEMPERATURE_LINE_NAME, ...
%                     "lineStyle",TEMPERATURE_LINE_STYLE, ...
%                     "lineWidth",TEMPERATURE_LINE_WIDTH, ...
%                     "labelName",TEMPERATURE_LABEL_NAME, ...
%                     "graphName",TEMPERATURE_GRAPH_TITLE, ...
%                     "lineWidth",[1,1,1,1], ..."yLimit",[22.4 27.5],...
%                     "isSave",FILE_IS_SAVE,"outDir",OUT_DIR_PATH, ...
%                     "fontSize",20,"LabelFontSize",30,"saveFileExt","jpg");
% %制御入力
% makeGraph(t',plotInputData , ...
%                     "lineName",CONTROLINPUT_LINE_NAME, ...
%                     "lineWidth",[1,1,1,1], ...
%                     "labelName",CONTROLINPUT_LABEL_NAME, ..."yLimit",[-0.8,1.6],...
%                     "location","southwest",...
%                     "graphName",CONTROLINPUT_GRAPH_TITLE, ...
%                     "isSave",FILE_IS_SAVE,"outDir",OUT_DIR_PATH, ...
%                     "fontSize",20,"LabelFontSize",30,"saveFileExt","jpg");













TEMPERATURE_GRAPH_TITLE = ["transmission time"];
TEMPERATURE_LINE_NAME = ["$T_{0}-r_1$","$T_{0}-r_3$","$\mathrm{Part} \mathrm{W_1}$","$$\mathrm{Part} \mathrm{W_3}$"];
TEMPERATURE_LINE_WIDTH = [1,1,2,2];
TEMPERATURE_LABEL_NAME = ["Time [$\mathrm{s}]$","TransTime [$\mathrm{s}]$"];
TEMPERATURE_LINE_STYLE = ["-"];


CONTROLINPUT_GRAPH_TITLE = ["execute time"];
CONTROLINPUT_LINE_NAME = ["$i_1$","$i_3$"];
CONTROLINPUT_LABEL_NAME = ["Time [$\mathrm{s}]$","ExecTime [$\mathrm{A}$]"];
%通信時間
makeGraph(t',timeDelayData, ..."lineName",TEMPERATURE_LINE_NAME, ...
                    "lineStyle",TEMPERATURE_LINE_STYLE, ..."lineWidth",TEMPERATURE_LINE_WIDTH, ...
                    "labelName",TEMPERATURE_LABEL_NAME, ...
                    "graphName",TEMPERATURE_GRAPH_TITLE, ...
                    "lineWidth",[1], ..."yLimit",[22.4 27.5],...
                    "isSave",FILE_IS_SAVE,"outDir",OUT_DIR_PATH, ...
                    "fontSize",20,"LabelFontSize",30,"saveFileExt","jpg");
%制御入力
makeGraph(t',execTimeData , ..."lineName",CONTROLINPUT_LINE_NAME, ...
                    "lineWidth",[1], ...
                    "labelName",CONTROLINPUT_LABEL_NAME, ..."yLimit",[-0.8,1.6],...
                    "location","southwest",...
                    "graphName",CONTROLINPUT_GRAPH_TITLE, ...
                    "isSave",FILE_IS_SAVE,"outDir",OUT_DIR_PATH, ...
                    "fontSize",20,"LabelFontSize",30,"saveFileExt","jpg");



