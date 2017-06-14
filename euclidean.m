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

