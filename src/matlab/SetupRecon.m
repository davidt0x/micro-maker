function RS = SetupRecon(params)
% SetupRecon - Function to create data structures needed to compute a
% 3D microstructure reconstruction from exemplar 2D data. This only sets
% up a reconstruction optimization for a single resolution. If you are
% looking for the multi-resolution analog (you most likely should be) then
% see SetupReconMultiRes.
%
% This function establishes data structures needed performing
% a 3D reconstruction optimization with specific parameters. This object
% can be passed to the Reconstruct function to perform this optimization.
%
% Syntax:  [RS] = SetupRecon(params)
%
% Inputs:
%    params - A struct specifying the reconstruction parameters
%
%       params.NB_SIZE: 
%           The neighborhood size to use. This is the size in pixels of the square
%           neighborhood. Smaller neighborhoods size should be used at
%           higher resolutions if reconstructions are too slow. NB_SIZE should be odd
%           and greater than 1.
%
%       params.RECON_SIZE: 
%           The dimensions of the reconstruction. Must be a postive vector
%           of length 3.
%
%       params.EXEMPLARS:
%           A cell array containing a single image (2D matrix of grayscale values
%           ranging between 0 and 1) for each element. These are the exemplar images
%           to use for the reconstructions. The orientations (within the
%           reconstruction) for these exemplars is defined by NB_INDICES. 
%
%       params.NB_INDICES:  
%           Currently, this code supports up to nine different exemplar image
%           orientations which are defined by the miller indices (0 0 1), (0 1 0),
%           (1 0 0), (1 1 0), (-1 1 0), (0 1 1), (0 -1 1), (1 0 1), and (-1 0 1). 
%           Consult, the code in MakeNBOffsets to understand more. NB_INDICES
%           defines, for each element of EXEMPLARS, which orientation the exemplar
%           has within reconstruction sample frame. Each element of NB_INDICES is an
%           index into the above list of nine orientations. 
%
%       params.NUM_CORES:  
%           The number of processor cores that the code can use.
%
%       params.ANN_ALGO
%           Either 'FLANN' or 'PatchTable'.
%
%       params.UseNeighborhoodSearchWeights
%           Whether to use the neighborhood weigthing scheme that tries to
%           promote diverse use of neughborhhods by downweighting overused
%           neighborhoods. By default we don't use this because it
%           increases search times (neighborhood search tables need to be
%           rebuilt) but it can help if you find the solution is getting
%           stuck in a local minimum by oversampling the same
%           neighborhoods.
%
% Outputs:
%    RS - The reconstruction struct containing all data structures
%    needed for performing a reconstruction optimization.
%
% Example: 
%    RS = SetupRecon(params);
%
% Other m-files required: 
% MAT-files required: MakeNBOffsets.m
%
% See also: SetupReconMultiRes, DestroyRecon

% Author: David Turner
% Email: davidt0x@gmail.com 
% Website: https://github.com/davidt0x
% December 2016

%------------- BEGIN CODE --------------

    % Add paths the the mex files if they aren't there. We will need it.
    add_paths;

    % Lets check an make sure we got a valid set of reconstruction parameters.
    [isValid, errMsg] = CheckReconParams(params);
    if(~isValid)
        error(errMsg);
    end

    % If the user didn't pass in the number of cores to use. Default to 1
    if(~isfield(params, 'NUM_CORES'))
        params.NUM_CORES = 1;
    end

    % If the user hasn't specified an ANN algorithm then use the default.
    if(~isfield(params, 'ANN_ALGO'))
       params.ANN_ALGO = 'FLANN'; 
    end

    % Copy all the parameters to the reconstruction object we will return.
    RS = params;

    % Compute half the neighborhood size. This is rounded down because we are
    % assuming the neighborhood is odd and we want to know how many pixels are
    % to the left and to the right of the center pixel.
    HALF_NB_SIZE = floor(RS.NB_SIZE/2);

    % Get the number of exemplars
    NUM_EXEMPLARS = length(RS.EXEMPLARS);

    % We need to pre-allocate some cell arays for the data structures.

    % This cell array contains, for each exemplar image, a matrix of all
    % extracted neighborhoods. That is, the matrix has "Number of Pixels"
    % rows and NB_SIZE*NB_SIZE columns.
    RS.Exemplar_NBs = cell(NUM_EXEMPLARS, 1);

    % Our reconstruction code uses approximate nearest neighbor lookup indices
    % from libraries like FLANN and PatchTable. Since these libraries are
    % implemented through MEX functions that create data structures in memory
    % outside of Matlab's environment we need to keep track of a pointer to
    % them that can be passed as a reference to these MEX calls. This pointer
    % is stored in this cell array, one for each exemplar image.
    RS.Exemplar_Index = cell(NUM_EXEMPLARS, 1);

    % For both FLANN and PatchTable there are important parameters to set that
    % determine their performance. These are the parameters set for each
    % exemplar image.
    RS.Exemplar_Params = cell(NUM_EXEMPLARS, 1);

    % Patches are extracted patches from each image. This lookup table
    % contains, for each pixel in an exemplar, the location of the extracted
    % patch centered on it in the Exemplar_NBs matrix.
    RS.Exemplar_NBLookup = cell(NUM_EXEMPLARS,1);

    % Patches are extracted patches from each image. This lookup table
    % contains, for each patch in Exemplar_NBs, the location of the pixel in
    % the exemplar that this patch is centred on.
    RS.NB_ExemplarLookup = cell(NUM_EXEMPLARS,1);

    % If the user passed in angle map for the images. Then store it
    isAngMap = isfield(RS, 'ANG_MAP');

    % If the exemplars passed in are non-binary images then we will
    % treat them as a grain map. We will store neighborhoods for
    % both the grain boundaries and the grain map itself.
    if(max(RS.EXEMPLARS{1}(:)) > 1)
        RS.Exemplar_NBs_Edges = cell(NUM_EXEMPLARS, 1);
        RS.E_Edges = RS.EXEMPLARS;
        for ii=1:NUM_EXEMPLARS
            RS.E_Edges{ii} = GetEdgesCircShift2D(RS.EXEMPLARS{ii});
        end
        RS.NB_GrainList = cell(NUM_EXEMPLARS, 1);
        RS.NB_NumGrains = cell(NUM_EXEMPLARS, 1);
    end

    % Get the neighborhoods offsets for this neighborhood size.
    % Use the standard list.
    nbOffsets = MakeNBOffsets(RS.NB_SIZE);
    nbOffsets = nbOffsets(RS.NB_INDICES);
    RS.nbOffsets = zeros([length(RS.NB_INDICES), size(nbOffsets{1})]);
    for ii=1:length(RS.NB_INDICES)
        RS.nbOffsets(ii, :) = nbOffsets{ii}(:);
    end

    % NNB Index Parameters for FLANN
    if(strcmp(RS.ANN_ALGO,'FLANN'))
        build_params.cores = RS.NUM_CORES;
        build_params.algorithm = 'kdtree';
        build_params.trees = 8;
        build_params.checks = 131;
    elseif(strcmp(RS.ANN_ALGO,'PatchTable'))
       build_params = [];
    else
        error('Invalid ANN_ALGO parameter. Must be either FLANN or PatchTable');
    end

    % Now, for each exemplar, lets construct the data structures.
    for ii=1:NUM_EXEMPLARS
        fprintf(1, '\tGetting neighborhoods for exemplar %d of %d ...', ii, NUM_EXEMPLARS);

        if(isfield(RS, 'NB_GrainList'))
            [NB_Hoods RS.NB_ExemplarLookup{ii} RS.NB_GrainList{ii}] = GetAllNHoodsMEX(RS.EXEMPLARS{ii}, RS.NB_SIZE);
            RS.NB_NumGrains{ii} = sum(RS.NB_GrainList{ii} > 0, 2);
            RS.NB_GrainList{ii} = uint16(RS.NB_GrainList{ii});
            RS.NB_GrainList{ii} = RS.NB_GrainList{ii}(:, 1:max(RS.NB_NumGrains{ii}));
        else
            [NB_Hoods RS.NB_ExemplarLookup{ii}] = GetAllNHoodsMEX(RS.EXEMPLARS{ii}, RS.NB_SIZE);
        end
        RS.Exemplar_NBs{ii} = NB_Hoods; 

        % If there is a grain boundary map. Then get the neighborhoods for this as well
        if(isfield(RS, 'E_Edges'))
            [RS.Exemplar_NBs_Edges{ii} RS.NB_ExemplarLookup{ii}] = GetAllNHoodsMEX(RS.E_Edges{ii}, RS.NB_SIZE);
        end

        % If we have an angle map. Then the images are indexed and we need to convert
        % the neighborhoods to their actual angles. We will use the 6,7,8 GSH coeffcients
        % because they have indpendent information.
        if(isAngMap)
            tmp = zeros(size(NB_Hoods, 1), 3*size(NB_Hoods, 2));
            for ll=1:size(ANG_MAP, 2)
                idx_start = (ll-1)*size(NB_Hoods, 2)+1;
                tmp(:, idx_start:(idx_start+size(NB_Hoods,2)-1)) = reshape(ANG_MAP(NB_Hoods, ll), size(NB_Hoods));
            end
            RS.Exemplar_NBs{ii} = tmp;
        end

        % Make a map for the original exemplar image that tells where each
        % neighborhood is stored in the the neighborhood table
        RS.Exemplar_NBLookup{ii} = zeros(size(RS.EXEMPLARS{ii}));
        for ll=1:size(RS.NB_ExemplarLookup{ii}, 1)
           RS.Exemplar_NBLookup{ii}(RS.NB_ExemplarLookup{ii}(ll, 1), RS.NB_ExemplarLookup{ii}(ll, 2)) = ll;  
        end

        fprintf(1, 'Done\n');

        % Build nearest neighbor query indices for these points.
        fprintf(1, '\tBuilding ANN index for exemplar %d of %d ...', ii, NUM_EXEMPLARS);
        RS.Exemplar_Params{ii} = build_params;
        tic;
        if(isfield(RS, 'Exemplar_NBs_Edges') && ~isAngMap)
            if(isfield(RS, 'NB_NumGrains'))
                if(strcmp(RS.ANN_ALGO, 'FLANN'))
                    [RS.Exemplar_Index{ii} RS.Exemplar_Params{ii}] = flann_build_index([RS.Exemplar_NBs_Edges{ii} RS.NB_NumGrains{ii}*1000]', build_params);
                elseif(strcmp(RS.ANN_ALGO, 'PatchTable'))
                    
                end
            else
                if(strcmp(RS.ANN_ALGO, 'FLANN'))
                    [RS.Exemplar_Index{ii} RS.Exemplar_Params{ii}] = flann_build_index(RS.Exemplar_NBs_Edges{ii}', build_params);
                elseif(strcmp(RS.ANN_ALGO, 'PatchTable'))
                    
                end
            end
        else
            if(strcmp(RS.ANN_ALGO, 'FLANN'))
                [RS.Exemplar_Index{ii} RS.Exemplar_Params{ii}] = flann_build_index([RS.Exemplar_NBs{ii}]' , build_params);
            elseif(strcmp(RS.ANN_ALGO, 'PatchTable'))
                
            end
        end

        elapsed_time = toc;
        fprintf(1, 'Done. Time to Build = %f seconds\n', elapsed_time);
    end

    % Setup the reconstruction change table. This tells us which
    % neighborhoods have changed in the reconstruction over the last
    % iteration. This allows us to only search for new matching
    % neighborhoods on those that have changed.
    RS.change_table = ones(RS.RECON_SIZE);

    % At each iteration we need to build a table that says for each voxel
    % in the reconstruction where the best matching neighborhoods is in each
    % exemplar.
    RS.NNB_Table = zeros([RS.RECON_SIZE NUM_EXEMPLARS]);

end