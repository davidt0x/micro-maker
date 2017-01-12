function P = BuildExemplarPyramid(E, num_levels)

if(nargin < 2)
    num_levels = 4;
end

P = cell(num_levels, 1);

P{1} = E;

for ii=2:num_levels
    P{ii} = downsample(P{ii-1});
    %P{ii} = double(round(P{ii}));
    %P{ii} = double(im2bw(P{ii}, graythresh(P{ii})));
end
    
