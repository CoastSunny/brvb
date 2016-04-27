function ico = renorm_forward_eeg(ico, o, varargin)

if (nargin > 2)
    for i=1:o
        ico(i).cfg = renorm_forward(ico(i).wh.mean, ico(i).wh.normal, 'eeg', varargin{1});
    end
else
    for i=1:o
        ico(i).cfg = renorm_forward(ico(i).wh.mean, ico(i).wh.normal, 'eeg');
    end
end
