%% plot_csv_cases.m
% 将不同 case 中相同指标画在同一张图上。
%
% 细线：t <= t_s
% 粗线：t >= t_s
%
% xlim:
%   [0, 所有 case 中最大 t_f]
%
% ylim:
%   根据每张图中所有 case 在 [t_s, t_f] 上的数据范围自动扩展 15%

clear; clc; close all;

%% ===== 用户可修改区 =====
dataDir = './csvs/crit';              % cases_data.mat 所在文件夹
matFile = fullfile(dataDir, "cases_data.mat");

t_s = 2;                    % 时间窗口起点，可自行修改，例如 t_s = 1000;

% 选择要画哪些图
doPlot.EHF = true;
doPlot.U = true;
doPlot.KE1_KE2 = true;
doPlot.PE = true;
doPlot.KE1_KE2_PE = false;
%% ======================

S = load(matFile, "Cases", "CaseTable");
Cases = S.Cases;

caseNames = fieldnames(Cases);

if isempty(caseNames)
    error("Cases 为空，请先运行 load_csv_cases.m。");
end

disp("将要绘制以下 cases：");
disp(caseNames);

if doPlot.EHF
    plot_metric_group(Cases, caseNames, ["EHF"], t_s, "EHF vs t");
end

if doPlot.U
    plot_metric_group(Cases, caseNames, ["U"], t_s, "U vs t");
end

if doPlot.KE1_KE2
    plot_metric_group(Cases, caseNames, ["KE1", "KE2"], t_s, "KE1 and KE2 vs t");
end

if doPlot.PE
    plot_metric_group(Cases, caseNames, ["PE"], t_s, "PE vs t");
end

if doPlot.KE1_KE2_PE
    plot_metric_group(Cases, caseNames, ["KE1", "KE2", "PE"], t_s, "KE1, KE2 and PE vs t");
end

%% ===== 局部函数 =====

function plot_metric_group(Cases, caseNames, yFields, t_s, figTitle)

    figure("Name", figTitle);
    hold on; grid on; box on;

    allYPost = [];
    maxTf = -inf;

    for i = 1:numel(caseNames)

        cname = caseNames{i};
        C = Cases.(cname);

        t = C.t(:);
        maxTf = max(maxTf, max(t));

        preMask = t <= t_s;
        postMask = t >= t_s;

        for j = 1:numel(yFields)

            yName = yFields(j);

            if ~isfield(C, yName)
                warning("case %s 中没有字段 %s，已跳过。", cname, yName);
                continue;
            end

            y = C.(yName);
            y = y(:);

            if numel(y) ~= numel(t)
                warning("case %s 中字段 %s 长度与 t 不一致，已跳过。", cname, yName);
                continue;
            end

            labelName = sprintf("%s.%s", cname, yName);

            % 先画粗线部分，保留 legend
            if any(postMask)
                h = plot(t(postMask), y(postMask), "-", ...
                    "LineWidth", 2.0, ...
                    "DisplayName", labelName);

                thisColor = h.Color;

                allYPost = [allYPost; y(postMask)];

                % 再画细线部分，不进入 legend
                if any(preMask)
                    plot(t(preMask), y(preMask), "-", ...
                        "LineWidth", 0.4, ...
                        "Color", thisColor, ...
                        "HandleVisibility", "off");
                end

            else
                % 如果 t_s 超过该 case 的最大时间，则整段都只画细线
                if any(preMask)
                    plot(t(preMask), y(preMask), "-", ...
                        "LineWidth", 0.4, ...
                        "DisplayName", labelName);
                end
            end
        end
    end

    xlabel("t");
    ylabel(strjoin(yFields, ", "));
    title(figTitle, "Interpreter", "none");

    if isfinite(maxTf) && maxTf > 0
        xlim([0, maxTf]);
    end

    % ylim 只根据 [t_s, t_f] 的粗线窗口确定
    if ~isempty(allYPost)
        yMin = min(allYPost, [], "omitnan");
        yMax = max(allYPost, [], "omitnan");

        if isfinite(yMin) && isfinite(yMax)
            if abs(yMax - yMin) < eps
                pad = max(0.1 * abs(yMin), 1e-12);
            else
                pad = 0.15 * (yMax - yMin);
            end

            ylim([yMin - pad, yMax + pad]);
        end
    end

    legend("Interpreter", "none", "Location", "best");
end
