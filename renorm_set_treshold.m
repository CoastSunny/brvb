function [nndips, ndip, histcut, hm, hv, treshold] = renorm_set_treshold(ndips, dip, varargin)

treshold = 1000;
if (nargin > 2 && varargin{1} > 0)
   
    histcut = varargin{1};
    
    [hm, hv] = hist(abs(dip), treshold);
    if (nargin > 3)
        ndip = find(abs(dip) <= hv(histcut)+hv(1));
    else
        ndip = find(abs(dip) >= hv(histcut)-hv(1));
    end
    
    if (histcut == 0)
        nndips = sum(hm);
    else
        nndips = sum(hm(histcut+1:end));
    end
    
else
    
    nndips = ndips+1;
    treshold = 1000;
    histcut = 2;
    [hm, hv] = hist(abs(dip), treshold);
    while nndips > ndips && histcut < 999
        ndip = find(abs(dip) >= hv(histcut));
        nndips = sum(hm(histcut:end));
        histcut = histcut + 1;
    end
    
end