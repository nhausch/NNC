function dmat = distance_matrix(Z)
[~, n] = size(Z);
dmat = zeros(n);
for i = 1:n
    for j = i:n
        if i ~= j
            dmat(j,i) = euclidean(Z(:,i), Z(:,j));
        end
    end
end
dmat = dmat+dmat';
end

