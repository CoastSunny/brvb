function ico = renorm_meandist(ico, O)

for o=1:O
    nfaces = size(ico(o).wh.faces, 1);
    
    dist=[];
    p=1;
    for i=1:nfaces
        x = ico(o).wh.faces(i,:);
        for j=1:size(x,2)
            [y ~] = find(ico(o).wh.faces == x(j));
            y(y==i) = [];
            for k = 1:size(y,1)
                dist(p) = sqrt(sum((ico(o).wh.mean(y(k),:) - ico(o).wh.mean(i,:)).^2)); %#ok<AGROW>
                p = p+1;
            end
        end
    end
    
    ico(o).wh.meandist = mean(dist);
    ico(o).wh.meandistvar = std(dist);
end