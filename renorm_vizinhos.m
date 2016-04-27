function [vizinhos1, vizinhos2] = renorm_vizinhos(faces)

nfaces = size(faces, 1);
vizinhos1 = zeros(nfaces, 3);
vizinhos2 = cell(nfaces, 1);
for i=1:nfaces
    x = faces(i,:);
    [a1 ~] = find(ismember(faces,x([1 2])));
    [a2 ~] = find(ismember(faces,x([1 3])));
    [a3 ~] = find(ismember(faces,x([2 3])));

    a4 = intersect(intersect(a1,a2),a3);
    a4(a4==i) = [];
    vizinhos1(i, :) = a4;

    a5 = unique([a1;a2;a3]);
    a5(a5==i) = [];
    a5(ismember(a5,a4)) = [];
    vizinhos2(i) = {a5};
end

