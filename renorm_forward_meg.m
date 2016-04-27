function ico = renorm_forward_meg(ico, o)

for i=1:o
    ico(i).cfg = renorm_forward(ico(i).wh.mean, ico(i).wh.normal, 'meg');
end