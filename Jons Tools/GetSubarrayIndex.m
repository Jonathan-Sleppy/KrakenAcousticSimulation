function [VertIndex,AdjacentVertIndex,HorizIndex,AdjacentHorizIndex] = GetSubarrayIndex(m,n,i)
    HorizIndex = (i-1)*m + 1:i*m;
    AdjacentHorizIndex = i*m +1:(i+1)*m;
    VertIndex = i:n:n*m;
    AdjacentVertIndex = i+1:n:n*m;
end
