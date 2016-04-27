function STATE = renorm_load2(datapath, varargin)

max = 1000;
if (nargin > 1)
    max = varargin{1};
end

alpha = false;
if (nargin > 2)
    alpha = varargin{2};
end

files = dir(datapath);
n = size(files,1);

if (n > max+2)
    n = max+2;
end

STATE = cell(1,n-2);
for i=1:n
    %    fprintf('opening %d...\n', i);
    fprintf('opening %d = %s...\n', i, files(i).name);
    [state , j] = renorm_load3(datapath, files(i), alpha);
    if ~isempty(state)
        STATE{j} = state;
    end
end
    