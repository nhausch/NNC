function [assignments,resolved,num_unassigned,numclusters] = min_neighbor_clusters(Z, k, num_neighbors)
[~, n] = size(Z);
maxiter = 100;
dmat = distance_matrix(Z);
assignments = zeros(1,n) - 1;
for i = 1:n
    [dvec,ivec] = sort(dmat(:,i));
    nproduct = 1;
    for j = 2:(num_neighbors+1)
        nproduct = nproduct * dvec(j);
    end
    sortvec(i) = nproduct;
    imat(:,i) = ivec;
end
imat(1,:) = [];
[~,order] = sort(sortvec);
numclusters = 0;
iter = 0;
resolved = false;
while ~resolved && iter < maxiter
    for i = 1:length(order)
        curr = order(i);
        nvec = imat(1:num_neighbors,curr);
        if assignments(curr) == -1
            punassigned = true;
        else
            punassigned = false;
            pindex = assignments(curr);
        end   
        cunassigned = true;
        for j = 1:num_neighbors
            if assignments(nvec(j)) ~= -1
                aindex = assignments(nvec(j));
                cunassigned = false;
            end
        end
        if numclusters < k && punassigned && cunassigned
            assignments(curr) = curr;
            for j = 1:num_neighbors
                assignments(nvec(j)) = curr;
            end
            numclusters = numclusters + 1;
        elseif punassigned && ~cunassigned
            for j = 1:num_neighbors
                assignments(nvec(j)) = aindex;
            end
            assignments(curr) = aindex;
        elseif ~punassigned && cunassigned
            for j = 1:num_neighbors
                assignments(nvec(j)) = assignments(curr);
            end
        elseif ~punassigned && ~cunassigned
            for j = 1:num_neighbors
                if assignments(nvec(j)) == -1
                    assignments(nvec(j)) = assignments(curr);
                else
                    if assignments(nvec(j)) ~= pindex
                        mergeindex = find(assignments == assignments(nvec(j)));
                        assignments(mergeindex) = pindex;
                        numclusters = numclusters - 1;
                    end
                end
            end
        end
    end
    if assignments(:) ~= -1
        resolved = true;
    end
iter = iter + 1;
end
temp = find(assignments == -1);
num_unassigned = length(temp);
xin = Z(1,:); yin = Z(2,:);
scatter(xin,yin,10,assignments);
end

function distance = euclidean(a,b)
sum = 0;
if length(a) == length(b)
    m = length(a);
else
    error('Vectors must have the same dimension.')
end
for i = 1:m
    sum = sum + abs(a(i)-b(i))^2;
end
distance = sqrt(sum);
end

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
