function [ a ] = find_first_common_ancestor( ancestors )
%FIND_LARGEST_COMMON_ANCESTOR Summary of this function goes here
%   every element is a list of ancestor

inter = ancestors{1};
for i=2:length(ancestors)
    inter = intersect( inter, ancestors{i} );
end

s = sort( inter );
a = s(1);
end

