function [ico, wh] = renorm_join(ico, rh, lh)

wh.orig.vertices = [lh.orig.vertices ; rh.orig.vertices];
wh.orig.faces = [lh.orig.faces ; rh.orig.faces + size(lh.orig.vertices,1)];

% wh.sphere.vertices = [lh.sphere.vertices ; rh.sphere.vertices];
d=100;
wh.sphere.vertices = [ [lh.sphere.vertices(:,1)-d lh.sphere.vertices(:,2) lh.sphere.vertices(:,3)]; [rh.sphere.vertices(:,1)+d rh.sphere.vertices(:,2) rh.sphere.vertices(:,3)] ];
wh.sphere.faces = [lh.sphere.faces ; rh.sphere.faces + size(lh.sphere.vertices,1)];

d=45;
wh.inflated.vertices = [ [lh.inflated.vertices(:,1)-d lh.inflated.vertices(:,2) lh.inflated.vertices(:,3)]; [rh.inflated.vertices(:,1)+d rh.inflated.vertices(:,2) rh.inflated.vertices(:,3)] ];
wh.inflated.faces = [lh.inflated.faces ; rh.inflated.faces + size(lh.inflated.vertices,1)];

wh.curv.cdata = [ lh.curv.cdata ;  rh.curv.cdata];

n = size(ico,2);

for i = 1:n

    ico(i).wh.mean = [ico(i).lh.mean ; ico(i).rh.mean];
    ico(i).wh.normal = [ico(i).lh.normal ; ico(i).rh.normal];
    ico(i).wh.faces = [ico(i).lh.faces ; ico(i).rh.faces + size(lh.orig.vertices,1)];
    ico(i).wh.aperface = [ico(i).lh.aperface ; ico(i).rh.aperface];
    ico(i).wh.vperface = [ico(i).lh.vperface ; ico(i).rh.vperface];

    if i > 1
        ico(i).wh.facesmap = [ico(i).lh.facesmap ; ico(i).rh.facesmap + size(ico(i-1).lh.faces,1)];
    else
        ico(1).wh.facesmap = [ico(1).lh.facesmap ; ico(1).rh.facesmap + size(ico(1).lh.faces,1)];
    end

    ico(i).wh.facesmapori = [ico(i).lh.facesmapori ; ico(i).rh.facesmapori + size(ico(i).lh.faces,1)];

    fprintf('building map to original surface %d ...  ', i);
    totalndips = length(ico(i).wh.faces);
    smap = cell(totalndips,1);
    for j=1:totalndips
        smap{j} = find(ico(i).wh.facesmapori==j);
    end
    ico(i).smap = smap;
    fprintf('done\n');
end

