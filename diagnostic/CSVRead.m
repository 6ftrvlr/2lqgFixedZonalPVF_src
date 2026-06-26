%% load_csv_cases.m
% 批量读取当前文件夹或指定文件夹中的 csv 文件，
% 自动判断 EHF 或 U 哪一列为常值，并将每个文件整理为一个 case 结构体。
%
% 输出：
%   Cases.(caseName)      : 结构体字典，推荐使用
%   CaseMap(caseName)     : containers.Map 字典，也可使用
%   CaseTable             : 汇总表
%
% 例：
%   load("cases_data.mat");
%   Cases.F_M0p25.t
%   Cases.F_M0p25.U
%   Cases.U_1p0.KE1

clear; clc;

%% ===== 用户可修改区 =====
dataDir = './csvs/samp';          % csv 文件所在文件夹；例如 "D:\data\csv_cases"
filePattern = "*.csv";

tol = 1e-15;            % 判断常值的误差限
nHead = 10;             % 先读前 10 行判断
nCheckMin = 100;        % 至少用后续 100 行校验

% 必须读入的字段。
% 如果以后要加字段，比如 "KE3"，就在这里追加：
% requiredVars = ["number", "t", "EHF", "U", "KE1", "KE2", "PE", "KE3"];
requiredVars = ["number", "t", "EHF", "U", "KE1", "KE2", "PE"];
%% ======================

files = dir(fullfile(dataDir, filePattern));

if isempty(files)
    error("没有在文件夹中找到 csv 文件：%s", dataDir);
end

Cases = struct();
CaseMap = containers.Map("KeyType", "char", "ValueType", "any");

CaseTable = table( ...
    strings(0,1), strings(0,1), zeros(0,1), strings(0,1), zeros(0,1), zeros(0,1), ...
    'VariableNames', ["name", "sourceFile", "type", "constantKind", "constantValue", "nRows"] ...
);

for k = 1:numel(files)

    filePath = fullfile(files(k).folder, files(k).name);
    fprintf("Reading %s\n", filePath);

    T = readtable(filePath, ...
        "FileType", "text", ...
        "VariableNamingRule", "preserve");

    varNames = string(T.Properties.VariableNames);

    missingVars = setdiff(requiredVars, varNames);
    if ~isempty(missingVars)
        error("文件 %s 缺少字段：%s", files(k).name, strjoin(missingVars, ", "));
    end

    nRows = height(T);

    if nRows < nHead + nCheckMin
        error("文件 %s 行数不足。至少需要 %d + %d 行。", ...
            files(k).name, nHead, nCheckMin);
    end

    C = struct();

    C.sourceFile = filePath;
    C.sourceName = files(k).name;
    C.nRows = nRows;

    % 读取所有要求字段为列向量
    for i = 1:numel(requiredVars)
        v = requiredVars(i);
        x = T.(v);

        if iscell(x) || isstring(x) || ischar(x)
            x = str2double(x);
        end

        C.(v) = x(:);
    end

    % 用前 10 行 + 后续至少 100 行校验 EHF / U 是否为常值
    idxHead = 1:nHead;
    idxCheck = (nHead + 1):nRows;   % 实际上使用了所有剩余行，强于只校验 100 行

    isEHFConst = is_constant_by_head_and_check(C.EHF, idxHead, idxCheck, tol);
    isUConst   = is_constant_by_head_and_check(C.U,   idxHead, idxCheck, tol);

    if isEHFConst && ~isUConst
        C.type = 1;
        C.constantKind = "EHF";
        C.constantValue = mean(C.EHF(idxHead), "omitnan");
        baseName = make_case_name("F", C.constantValue);

    elseif isUConst && ~isEHFConst
        C.type = 0;
        C.constantKind = "U";
        C.constantValue = mean(C.U(idxHead), "omitnan");
        baseName = make_case_name("U", C.constantValue);

    elseif isEHFConst && isUConst
        error("文件 %s 中 EHF 和 U 都被判定为常值，命名存在歧义。", files(k).name);

    else
        error("文件 %s 中 EHF 和 U 都没有被判定为常值。", files(k).name);
    end

    % 若出现重复 case 名，自动追加编号，避免覆盖
    caseName = baseName;
    rep = 2;
    while isfield(Cases, caseName)
        caseName = sprintf("%s_rep%d", baseName, rep);
        rep = rep + 1;
    end

    C.name = caseName;
    C.baseName = baseName;

    % 一些方便后处理的附加字段
    C.tStart = C.t(1);
    C.tEnd = C.t(end);

    Cases.(caseName) = C;
    CaseMap(char(caseName)) = C;

    CaseTable = [CaseTable; {
        string(caseName), string(files(k).name), C.type, C.constantKind, C.constantValue, C.nRows
    }];
end

save(fullfile(dataDir, "cases_data.mat"), "Cases", "CaseMap", "CaseTable");

disp("读取完成。汇总如下：");
disp(CaseTable);

%% ===== 局部函数 =====

function tf = is_constant_by_head_and_check(x, idxHead, idxCheck, tol)
    x = x(:);

    x0 = mean(x(idxHead), "omitnan");

    headOK = max(abs(x(idxHead) - x0), [], "omitnan") <= tol;
    checkOK = max(abs(x(idxCheck) - x0), [], "omitnan") <= tol;

    tf = headOK && checkOK;
end

function name = make_case_name(prefix, value)
    % 将数值转换为合法 MATLAB 字段名。
    %
    % 例：
    %   EHF = -0.25  -> F_M0p25
    %   U   = 1.0    -> U_1
    %   U   = -2.5   -> U_M2p5
    %
    % 小数点用 p 表示，负号用 M 表示。

    if value < 0
        signPart = "M";
    else
        signPart = "";
    end

    valAbs = abs(value);

    numStr = sprintf("%.15g", valAbs);
    numStr = strrep(numStr, ".", "p");
    numStr = strrep(numStr, "e-", "em");
    numStr = strrep(numStr, "e+", "ep");
    numStr = strrep(numStr, "+", "");

    rawName = prefix + "_" + signPart + numStr;

    % 双保险：保证一定是合法字段名
    name = matlab.lang.makeValidName(rawName);
end
