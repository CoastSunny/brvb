function STATE = renorm_load(datapath, varargin)

global DATAPATH;
global SUBJ;
global ANALISES;

if nargin < 1
    datapath = [DATAPATH filesep SUBJ filesep ANALISES ];
end

max = 50;
if (nargin > 1)
    max = varargin{1};
end

files = dir(datapath);
n = size(files,1);

if (n > max+2)
    n = max+2;
end

STATE = cell(1,n-2);
j=1;
for i=1:n
    %    fprintf('opening %d...\n', i);
%     if str2num(files(i).name)
    if files(i).isdir && ~strcmp(files(i).name, '.') && ~strcmp(files(i).name, '..')
        fprintf('opening %d = %s...\n', i, files(i).name);
        currfile = [datapath filesep files(i).name filesep 'STATE.mat'];
%         STATE(str2num(files(i).name)) = {open(path)};
%         STATE{str2num(files(i).name)}.filename = files(i).name;
        STATE(j) = {open(currfile)};
        STATE{j}.filename = files(i).name;
        j=j+1;
    end
end
    
