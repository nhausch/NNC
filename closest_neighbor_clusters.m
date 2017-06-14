function [assignments, numclusters, num_neighbors] = closest_neighbor_clusters(Z, k, sigma)
    dmat = distance_matrix(Z);
    [~, n] = size(Z);
    for i = 1:n
        [dvec,ivec] = sort(dmat(:,i));
        sdmat(:,i) = dvec;
        imat(:,i) = ivec;
    end
    for num_neighbors = 3:100
        [assignments,resolved,~,numclusters] = min_neighbor_clusters_2(Z, dmat, imat, k, num_neighbors, sigma);
        if resolved == 1
            break;
        end
    end
end

function [assignments,resolved,num_unassigned,numclusters] = min_neighbor_clusters_2(Z, sdmat, imat, k, num_neighbors, sigma)
[~, n] = size(Z);
thresh = floor(n*sigma);
maxiter = 100;
assignments = zeros(1,n) - 1;
for i = 1:n
    nproduct = 1;
    for j = 2:(num_neighbors+1)
        nproduct = nproduct * sdmat(j,i);
    end
    sortvec(i) = nproduct;
end
imat(1,:) = [];
[~,order] = sort(sortvec);
numclusters = 0;
iter = 0;
resolved = false;
while ~resolved && iter < maxiter
    for i = 1:length(order)
        if i > thresh
            isout = true;
        else
            isout = false;
        end
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
        if numclusters < k && punassigned && cunassigned && ~isout
            assignments(curr) = curr;
            for j = 1:num_neighbors
                assignments(nvec(j)) = curr;
            end
            numclusters = numclusters + 1;
        elseif punassigned && ~cunassigned
            assignments(curr) = aindex;
            if ~isout
                for j = 1:num_neighbors
                    assignments(nvec(j)) = aindex;
                end
            end
        elseif ~punassigned && cunassigned
            for j = 1:num_neighbors
                assignments(nvec(j)) = assignments(curr);
            end
        elseif ~punassigned && ~cunassigned
            for j = 1:num_neighbors
                if assignments(nvec(j)) == -1
                    assignments(nvec(j)) = assignments(curr);
                else
                    if ~isout
                        if assignments(nvec(j)) ~= pindex
                            mergeindex = find(assignments == assignments(nvec(j)));
                            assignments(mergeindex) = pindex;
                            numclusters = numclusters - 1;
                        end
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
