function [nbOffsets] = MakeNBOffsets(NB_SIZE)
% This function makes a standard set of 3D neighborhood offsets
% It produces 9 at the moment.

    HALF_NB_SIZE = floor(NB_SIZE/2);
    HS = HALF_NB_SIZE;

    % Generate the neighborhoods offsets for the three orthogonal
    % neighborhoods that are aligned with rows, columns, and layers of the 
    % 3D image.
    [x y z] = ndgrid((0-HS):(0+HS), (0-HS):(0+HS), 0); 
    nbOffsets{1} = [x(:) y(:) z(:)]; % (0,0,1)
    nbOffsets{2} = [z(:) x(:) y(:)]; % (0,1,0)
    nbOffsets{3} = [x(:) z(:) y(:)]; % (1,0,0)

    % Now get the offsets for the diagonal neighborhoods that are 45
    % degrees between any two of the previous neighborhoods
    w = repmat( ((0-HS):(0+HS))', [NB_SIZE, 1]);
    nbOffsets{4} = [w flipud(x(:)) y(:)];       % (1,1,0)
    nbOffsets{5} = [w x(:) y(:)];               % (-1,1,0)
    nbOffsets{6} = [x(:) y(:) flipud(w)];       % (0,1,1)
    nbOffsets{7} = [x(:) y(:) w];               % (0,-1,1)
    nbOffsets{8} = [x(:) y(:) flipud(y(:))];    % (1,0,1)
    nbOffsets{9} = [x(:) y(:) y(:)];            % (-1,0,1)
    
    %for ii=1:3
    %     M = zeros(10,10,10); 
    %     T = nbOffsets{ii}; 
    %     M(sub2ind(size(M), T(:,1)+5, T(:,2)+5, T(:,3)+5)) = ii; 
    %     isosurface(M, ii-1);
    %     axis equal; axis tight; box on; grid on;
    %     xlabel('X');
    %     ylabel('Y');
    %     zlabel('Z');
    %end
end
