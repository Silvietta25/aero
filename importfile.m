function polareNaca0012 = importfile(filename, dataLines)
%IMPORTFILE Import data from a text file
%  POLARENACA0012 = IMPORTFILE(FILENAME) reads data from text file
%  FILENAME for the default selection.  Returns the numeric data.
%
%  POLARENACA0012 = IMPORTFILE(FILE, DATALINES) reads data for the
%  specified row interval(s) of text file FILENAME. Specify DATALINES as
%  a positive scalar integer or a N-by-2 array of positive scalar
%  integers for dis-contiguous row intervals.
%
%  Example:
%  polareNaca0012 = importfile("/home/luigi/Hess Smith method/Hess Smith classico/polareNaca0012.txt", [13, Inf]);
%
%  See also READTABLE.
%
% Auto-generated by MATLAB on 19-Dec-2023 18:42:49

%% Input handling

% If dataLines is not specified, define defaults
if nargin < 2
    dataLines = [13, Inf];
end

%% Set up the Import Options and import the data
opts = delimitedTextImportOptions("NumVariables", 13);

% Specify range and delimiter
opts.DataLines = dataLines;
opts.Delimiter = " ";

% Specify column names and types
opts.VariableNames = ["VarName1", "VarName2", "Var3", "Var4", "Var5", "Var6", "Var7", "Var8", "Var9", "Var10", "Var11", "Var12", "Var13"];
opts.SelectedVariableNames = ["VarName1", "VarName2"];
opts.VariableTypes = ["double", "double", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string"];

% Specify file level properties
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";
opts.ConsecutiveDelimitersRule = "join";
opts.LeadingDelimitersRule = "ignore";

% Specify variable properties
opts = setvaropts(opts, ["Var3", "Var4", "Var5", "Var6", "Var7", "Var8", "Var9", "Var10", "Var11", "Var12", "Var13"], "WhitespaceRule", "preserve");
opts = setvaropts(opts, ["Var3", "Var4", "Var5", "Var6", "Var7", "Var8", "Var9", "Var10", "Var11", "Var12", "Var13"], "EmptyFieldRule", "auto");
opts = setvaropts(opts, "VarName1", "TrimNonNumeric", true);
opts = setvaropts(opts, "VarName1", "ThousandsSeparator", ",");

% Import the data
polareNaca0012 = readtable(filename, opts);

%% Convert to output type
polareNaca0012 = table2array(polareNaca0012);
end