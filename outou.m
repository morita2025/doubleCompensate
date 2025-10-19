%-------------------------------------------------------------
% パラメータ1次遅れシステムのステップ応答比較プログラム
%-------------------------------------------------------------
clear; clc; close all;

% ゲイン（固定:1）
K = 1;

% 時定数の候補（自由に追加可能）
tau_list = [1/10, 1/5, 1/2, 1];   % [s]

% 図を準備
figure;
hold on;
grid on;

% それぞれの時定数に対して伝達関数を定義し、ステップ応答を描画
for tau = tau_list
    % 伝達関数 G(s) = K / (tau*s + 1)
    sys = tf(K, [tau 1]);
    
    % ステップ応答を描画
    step(sys);
end

% 凡例を追加
legendStrings = arrayfun(@(x) sprintf('\\tau = %.1f', x), tau_list, 'UniformOutput', false);
legend(legendStrings, 'Location', 'southeast');

% 軸ラベル・タイトル
xlabel('時間 [s]');
ylabel('出力');
title('1次遅れシステムのステップ応答比較');
