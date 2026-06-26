%% CasStat.m
% 根据若干诊断量序列判断统计稳态开始时间，并计算稳态段统计量
% 输出结果追加保存到 cases_data.mat

clear; clc;

%% ===================== 用户设置区 =====================
dataDir = './csvs/crit';              % cases_data.mat 所在文件夹
matFile = fullfile(dataDir, "cases_data.mat");

% 待检测的诊断量，也就是 struct 中标注为 o 的序列
p.diagFields = ["U"]; % Lite
%p.diagFields = ["EHF", "U", "KE1", "KE2", "PE"]; % Too ideal.

% 时间窗判据参数
p.windowLength = 5;          % 时间窗长度 W
p.windowStep   = 0.02;         % 窗口移动步长
p.meanTol      = 0.10;         % 均值偏移阈值，例如 1 %
p.sdTol        = 0.10;         % SD 偏移阈值，例如 3 %
p.nPass        = 10;          % 连续通过次数
p.tEarliest    = 3.0;          % 最早允许判断稳定时间，例如 t = 4

% 归一化尺度选项：
% "postTMean"：均值尺度取 t >= p.refTime 后全体值的均值绝对值；
%              SD 尺度取 t >= p.refTime 后全体值的标准差
% "none"     ：不设参考值，均值尺度取当前后窗均值绝对值；
%              若均值接近零，程序报错
p.refMode = "postTMean";
p.refTime = 3.0;
p.zeroTol = 1e-12;

% 作图选择："meanStd", "box", "both", "none"
p.plotMode = "both";

% pointType 字段：先设两个，后面 assignPointTypes() 中可自行修改
p.pointTypeFields = ["pointType1", "pointType2"];

% 是否保存图片
p.saveFigures = false;
p.figDir = fullfile(dataDir, "stat_figures");

%% ===================== 读入数据 =====================
S = load(matFile);

if ~isfield(S, 'CaseMap')
    error('cases_data.mat 中未找到 CaseMap。');
end

CaseMap = S.CaseMap;

if isfield(S, 'cns')
    cns = string(S.cns(:));
else
    cns = string(keys(CaseMap));
end

CaseStatMap = containers.Map('KeyType', 'char', 'ValueType', 'any');

Rows = cell(0, 21);
rowNames = { ...
    'cn', 'field', 'constantAbs', 'constantKind', 'caseName', 'baseName', ...
    'pointType1', 'pointType2', ...
    'TstSystem', 'TstField', ...
    'nSample', 'meanValue', 'sdValue', 'varianceValue', ...
    'skewnessValue', 'kurtosisValue', ...
    'minValue', 'q25', 'q50', 'q75', 'maxValue'};

%% ===================== 主循环 =====================
for ic = 1:numel(cns)
    cn = char(cns(ic));
    C = CaseMap(cn);

    if ~isfield(C, 't')
        error('case %s 中没有字段 t。', cn);
    end

    t = double(C.t(:));
    if any(diff(t) <= 0)
        error('case %s 的 t 不是严格递增序列。', cn);
    end

    [pointType1, pointType2] = assignPointTypes(C);

    TstByField = struct();
    TestDetail = struct();

    fprintf('\nCase %s\n', cn);

    % -------- 对每个诊断量分别判稳 --------
    for jf = 1:numel(p.diagFields)
        field = char(p.diagFields(jf));

        if ~isfield(C, field)
            error('case %s 中没有诊断量字段 %s。', cn, field);
        end

        X = double(C.(field)(:));
        if numel(X) ~= numel(t)
            error('case %s 中 %s 与 t 长度不一致。', cn, field);
        end

        [TstField, detail] = findSteadyTimeByWindows(t, X, p, cn, field);

        TstByField.(field) = TstField;
        TestDetail.(field) = detail;

        fprintf('  %-6s T_st = %.8g\n', field, TstField);
    end

    % -------- 系统稳态时间：所有诊断量判稳时间的最大值 --------
    TstValues = struct2array(TstByField);
    TstSystem = max(TstValues);

    fprintf('  System T_st = %.8g\n', TstSystem);

    % -------- 统计稳态段统计量 --------
    CaseStat = struct();
    CaseStat.cn = string(cn);
    CaseStat.constantAbs = abs(getNumericScalar(C, 'constantValue'));
    CaseStat.constantKind = getFieldAsString(C, 'constantKind');
    CaseStat.caseName = getFieldAsString(C, 'name');
    CaseStat.baseName = getFieldAsString(C, 'baseName');
    CaseStat.pointType1 = pointType1;
    CaseStat.pointType2 = pointType2;
    CaseStat.TstSystem = TstSystem;
    CaseStat.TstByField = TstByField;
    CaseStat.TestDetail = TestDetail;
    CaseStat.stats = struct();

    for jf = 1:numel(p.diagFields)
        field = char(p.diagFields(jf));
        X = double(C.(field)(:));

        maskSteady = t >= TstSystem;
        Xs = X(maskSteady);

        st = calcStatsFirstFourAndBox(Xs);
        CaseStat.stats.(field) = st;

        Rows(end+1, :) = { ...
            string(cn), string(field), CaseStat.constantAbs, ...
            CaseStat.constantKind, CaseStat.caseName, CaseStat.baseName, ...
            pointType1, pointType2, ...
            TstSystem, TstByField.(field), ...
            st.n, st.meanValue, st.sdValue, st.varianceValue, ...
            st.skewnessValue, st.kurtosisValue, ...
            st.minValue, st.q25, st.q50, st.q75, st.maxValue}; %#ok<SAGROW>
    end

    % 把系统稳态时间也写回原 CaseMap，便于后续调用
    C.T_st = TstSystem;
    C.T_st_byField = TstByField;
    CaseMap(cn) = C;

    CaseStatMap(cn) = CaseStat;
end

%% ===================== 整理为表格并保存 =====================
SteadyStatTable = cell2table(Rows, 'VariableNames', rowNames);

strCols = ["cn", "field", "constantKind", "caseName", "baseName", ...
           "pointType1", "pointType2"];
for k = 1:numel(strCols)
    SteadyStatTable.(strCols(k)) = string(SteadyStatTable.(strCols(k)));
end

save(matFile, 'CaseMap', 'CaseStatMap', 'SteadyStatTable', 'p', '-append');

fprintf('\n统计稳态后处理已保存到：\n%s\n', matFile);

%% ===================== 作图 =====================
makePlots(SteadyStatTable, CaseMap, CaseStatMap, p);

%% ========================================================================
%%                               局部函数
%% ========================================================================

function [Tst, detail] = findSteadyTimeByWindows(t, X, p, cn, field)
    W = p.windowLength;
    dW = p.windowStep;

    if t(end) - t(1) < 2*W
        error('case %s, field %s: 时间长度不足两个窗口。', cn, field);
    end

    ref = getReferenceScales(t, X, p, cn, field);

    % 对相邻窗口 W1=[s,s+W), W2=[s+W,s+2W) 做比较；
    % 稳态时间定义为最后一个通过判断的 W2 的前端，即 s+W。
    sFirst = max(t(1), p.tEarliest - W);
    sLast  = t(end) - 2*W;

    starts = sFirst:dW:sLast;

    if isempty(starts)
        error('case %s, field %s: 可用于判稳的窗口为空。', cn, field);
    end

    rec = nan(numel(starts), 9);
    runPass = 0;
    Tst = NaN;

    for k = 1:numel(starts)
        s = starts(k);

        mask1 = (t >= s)     & (t < s + W);
        mask2 = (t >= s + W) & (t < s + 2*W);

        x1 = X(mask1);
        x2 = X(mask2);

        x1 = x1(isfinite(x1));
        x2 = x2(isfinite(x2));

        if numel(x1) < 2 || numel(x2) < 2
            passNow = false;
            meanRel = NaN;
            sdRel = NaN;
            mu1 = NaN; mu2 = NaN;
            sd1 = NaN; sd2 = NaN;
        else
            mu1 = mean(x1);
            mu2 = mean(x2);
            sd1 = std(x1, 0);
            sd2 = std(x2, 0);

            if strcmpi(char(p.refMode), 'postTMean')
                meanScale = ref.meanScale;
                sdScale   = ref.sdScale;
            elseif strcmpi(char(p.refMode), 'none')
                meanScale = abs(mu2);
                sdScale   = sd2;
            else
                error('未知 p.refMode：%s', p.refMode);
            end

            if meanScale < p.zeroTol
                error(['case %s, field %s: 均值归一化尺度接近零。', ...
                       '请改用参考尺度，或确认该变量不适合用相对均值偏移判据。'], ...
                       cn, field);
            end

            meanRel = abs(mu2 - mu1) / meanScale;

            sdDiff = abs(sd2 - sd1);
            if sdScale < p.zeroTol
                if sdDiff < p.zeroTol
                    sdRel = 0;
                else
                    error('case %s, field %s: SD 归一化尺度接近零但 SD 偏移非零。', cn, field);
                end
            else
                sdRel = sdDiff / sdScale;
            end

            passNow = (meanRel <= p.meanTol) && (sdRel <= p.sdTol);
        end

        if passNow
            runPass = runPass + 1;
        else
            runPass = 0;
        end

        front2 = s + W;

        rec(k, :) = [front2, mu1, mu2, sd1, sd2, meanRel, sdRel, double(passNow), runPass];

        if runPass >= p.nPass
            Tst = front2;
            break;
        end
    end

    if isnan(Tst)
        error('case %s, field %s: 未找到满足连续通过次数的统计稳态时间。', cn, field);
    end

    rec = rec(1:k, :);
    detail = struct();
    detail.Tst = Tst;
    detail.windowLength = W;
    detail.windowStep = dW;
    detail.meanTol = p.meanTol;
    detail.sdTol = p.sdTol;
    detail.nPass = p.nPass;
    detail.passTable = array2table(rec, 'VariableNames', ...
        {'Tfront2','meanWin1','meanWin2','sdWin1','sdWin2', ...
         'meanRelShift','sdRelShift','passNow','runPass'});
end

function ref = getReferenceScales(t, X, p, cn, field)
    if strcmpi(char(p.refMode), 'postTMean')
        mask = t >= p.refTime;
        xref = X(mask);
        xref = xref(isfinite(xref));

        if numel(xref) < 2
            error('case %s, field %s: t >= %.8g 后样本数不足。', cn, field, p.refTime);
        end

        ref.meanScale = abs(mean(xref));
        ref.sdScale = std(xref, 0);

        if ref.meanScale < p.zeroTol
            error('case %s, field %s: t >= %.8g 后均值接近零。', cn, field, p.refTime);
        end
    else
        ref.meanScale = NaN;
        ref.sdScale = NaN;
    end
end

function st = calcStatsFirstFourAndBox(X)
    X = double(X(:));
    X = X(isfinite(X));

    if isempty(X)
        error('稳态段样本为空，无法计算统计量。');
    end

    mu = mean(X);
    sd = std(X, 0);
    varianceValue = sd^2;

    if sd <= 0
        skewnessValue = NaN;
        kurtosisValue = NaN;
    else
        z = (X - mu) / sd;
        skewnessValue = mean(z.^3);
        kurtosisValue = mean(z.^4);
    end

    q = localQuantile(X, [0.25, 0.50, 0.75]);

    st = struct();
    st.n = numel(X);
    st.meanValue = mu;
    st.sdValue = sd;
    st.varianceValue = varianceValue;
    st.skewnessValue = skewnessValue;
    st.kurtosisValue = kurtosisValue;
    st.minValue = min(X);
    st.q25 = q(1);
    st.q50 = q(2);
    st.q75 = q(3);
    st.maxValue = max(X);
end

function q = localQuantile(X, probs)
    X = sort(X(:));
    n = numel(X);
    q = nan(size(probs));

    if n == 1
        q(:) = X(1);
        return;
    end

    for i = 1:numel(probs)
        p0 = probs(i);
        pos = 1 + (n - 1) * p0;
        lo = floor(pos);
        hi = ceil(pos);

        if lo == hi
            q(i) = X(lo);
        else
            w = pos - lo;
            q(i) = (1 - w) * X(lo) + w * X(hi);
        end
    end
end

function [pointType1, pointType2] = assignPointTypes(C)
    % ===== 用户可修改区：给不同 case 打标签，方便后续按类别作图 =====
    % 当前默认：
    % pointType1 = type
    % pointType2 = constantKind

    pointType1 = getFieldAsString(C, 'type');
    pointType2 = getFieldAsString(C, 'constantKind');

    if strlength(pointType1) == 0
        pointType1 = "undefinedType";
    end

    if strlength(pointType2) == 0
        pointType2 = "undefinedConstantKind";
    end

    % 预留添加字段位置：
    % 例如：
    % pointType3 = ...
    % pointType4 = ...
end

function val = getNumericScalar(C, fieldName)
    if ~isfield(C, fieldName)
        val = NaN;
        return;
    end

    v = C.(fieldName);

    if isnumeric(v)
        val = double(v(1));
    elseif isstring(v) || ischar(v)
        val = str2double(string(v));
    else
        val = NaN;
    end
end

function s = getFieldAsString(C, fieldName)
    if ~isfield(C, fieldName)
        s = "";
        return;
    end

    v = C.(fieldName);

    if isstring(v)
        s = v(1);
    elseif ischar(v)
        s = string(v);
    elseif isnumeric(v) && isscalar(v)
        s = string(v);
    else
        s = string(class(v));
    end
end

function makePlots(StatTable, CaseMap, CaseStatMap, p)
    mode = lower(char(p.plotMode));

    if strcmp(mode, 'none')
        return;
    end

    if p.saveFigures && ~exist(p.figDir, 'dir')
        mkdir(p.figDir);
    end

    if strcmp(mode, 'meanstd') || strcmp(mode, 'both')
        plotMeanStd(StatTable, p);
    end

    if strcmp(mode, 'box') || strcmp(mode, 'both')
        plotBoxByConstantValue(StatTable, CaseMap, CaseStatMap, p);
    end
end

function plotMeanStd(StatTable, p)
    for jf = 1:numel(p.diagFields)
        field = p.diagFields(jf);
        T = StatTable(StatTable.field == field, :);

        if isempty(T)
            continue;
        end

        [~, ord] = sort(T.constantAbs);
        T = T(ord, :);

        figure('Name', sprintf('%s mean ± SD', field));
        hold on;

        types = unique(T.pointType1, 'stable');

        for it = 1:numel(types)
            idx = T.pointType1 == types(it);

            errorbar(T.constantAbs(idx), ...
                     T.meanValue(idx), ...
                     T.sdValue(idx), ...
                     'o-', ...
                     'DisplayName', char(types(it)));
        end

        xlabel('|constantValue|');
        ylabel(sprintf('%s: mean ± SD', field));
        title(sprintf('%s versus |constantValue|', field));
        grid on;

        if numel(types) > 1
            legend('Location', 'best');
        end

        if p.saveFigures
            saveas(gcf, fullfile(p.figDir, sprintf('%s_meanSD.png', char(field))));
        end
    end
end

function plotBoxByConstantValue(StatTable, CaseMap, CaseStatMap, p)
    for jf = 1:numel(p.diagFields)
        field = p.diagFields(jf);
        T = StatTable(StatTable.field == field, :);

        if isempty(T)
            continue;
        end

        xUnique = unique(T.constantAbs);
        values = [];
        groups = [];

        for r = 1:height(T)
            cn = char(T.cn(r));
            C = CaseMap(cn);
            CS = CaseStatMap(cn);

            t = double(C.t(:));
            X = double(C.(char(field))(:));

            mask = t >= CS.TstSystem;
            xs = X(mask);
            xs = xs(isfinite(xs));

            g = find(xUnique == T.constantAbs(r), 1, 'first');

            values = [values; xs(:)]; %#ok<AGROW>
            groups = [groups; g * ones(numel(xs), 1)]; %#ok<AGROW>
        end

        figure('Name', sprintf('%s box plot', field));

        if exist('boxchart', 'file') == 2
            boxchart(groups, values);
            xticks(1:numel(xUnique));
            xticklabels(compose('%.4g', xUnique));
        else
            boxplot(values, groups, 'Labels', cellstr(compose('%.4g', xUnique)));
        end

        xlabel('|constantValue|');
        ylabel(char(field));
        title(sprintf('%s steady-state samples by |constantValue|', field));
        grid on;

        if p.saveFigures
            saveas(gcf, fullfile(p.figDir, sprintf('%s_box.png', char(field))));
        end
    end
end
