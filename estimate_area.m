function A = estimate_area(wh, ico)

A = zeros(1,6);

finer = size(ico,2);

for o=1:finer
    n = size(ico(o).wh.faces, 1);
    A(o) = 0;
    max=10;
    for i = 1:n
        v13 = wh.orig.vertices(ico(o).wh.faces(i,1), :) - wh.orig.vertices(ico(o).wh.faces(i,3), :);
        v12 = wh.orig.vertices(ico(o).wh.faces(i,1), :) - wh.orig.vertices(ico(o).wh.faces(i,2), :);
        A(o) = A(o) + sqrt(sum(cross(v12, v13).^2))/2;
        if (rem(i,max) == 0)
            fprintf('new : %d de %d\n', i, n);
        end
    end
    fprintf('Area para ordem %d : %f cm^2 (media de %f mm^2 por face) \n', o, A(o)/100, A(o)/n);
end

n = size(wh.orig.faces, 1);
A(6) = 0;
max=ceil(n/1000);
for i = 1:n
    v13 = wh.orig.vertices(wh.orig.faces(i,1), :) - wh.orig.vertices(wh.orig.faces(i,3), :);
    v12 = wh.orig.vertices(wh.orig.faces(i,1), :) - wh.orig.vertices(wh.orig.faces(i,2), :);
%    A(6) = A(6) + sum(sqrt(cross(v12, v13).^2))/2;
    A(6) = A(6) + sqrt(sum(cross(v12, v13).^2))/2;
    if (rem(i,max) == 0)
        fprintf('%d de %d\n', i, n);
    end
end
fprintf('Area para supreficie original : %f cm^2 (media de %f mm^2 por face)\n', A(6)/100, A(6)/n);

