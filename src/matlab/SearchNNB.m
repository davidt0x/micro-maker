function [Recon] = SearchNNB(S, Recon, useWeights)

    if(nargin < 3)
        useWeights = 0;
    end

    if(max(S(:)) >= 2 && ~isfield(Recon, 'ANG_MAP'))
        Sgrain = S;
        S = GetEdgesCircShift(S);
    end

    % For each exemplar, we need to extract neighborhoods from the 3D image
    % along the same orientation and find the best matching neighborhoods
    % in the 2D exemplar.
    for ExIndex=1:size(Recon.EXEMPLARS, 1)
        nbOffsets = squeeze(Recon.nbOffsets(ExIndex, :, :));
        [NB_queries, isPeriodic] = Extract3DNeighborhoods(S, nbOffsets, Recon.NUM_CORES);
        
        % If there is an angle map then we need to convert the indexed image to
        % actual values
        if(isfield(Recon, 'ANG_MAP'))
            ANG_MAP = Recon.ANG_MAP;
            tmp = zeros(size(NB_queries, 1), 3*size(NB_queries, 2));
            for ii=1:size(ANG_MAP, 2)
                idx_start = (ii-1)*size(NB_queries, 2)+1;
                tmp(:, idx_start:(idx_start+size(NB_queries,2)-1)) = reshape(ANG_MAP(NB_queries, ii), size(NB_queries));
            end
            NB_queries = tmp;
        end

        % If we are using neighborhood weights, then we need to append them
        % to the neighborhoods themselves before doing our query. They will
        % just be an additional dimension.
        if(useWeights)
            tmp = zeros(size(NB_queries,1), size(NB_queries,2)+1);
            tmp(:, 1:end-1) = NB_queries;
            NB_queries = tmp;
            
            NB_db = zeros(size(Recon.Exemplar_NBs{ExIndex},1), size(Recon.Exemplar_NBs{ExIndex},2)+1);
            if(isfield(Recon, 'Exemplar_NBs_Edges'))
                NB_db(:, 1:end-1) = Recon.Exemplar_NBs_Edges{ExIndex};
            else
                NB_db(:, 1:end-1) = Recon.Exemplar_NBs{ExIndex};
            end
                
            NB_db(:, end) = Recon.NBWeights{ExIndex};
                
            [Xnnidx] = flann_search(NB_db', NB_queries', 1, Recon.Exemplar_Params{ExIndex});
        else
            if(isfield(Recon, 'NB_NumGrains') && ~isfield(Recon, 'ANG_MAP'))
                [NB_queries_grain isPeriodic GrainList] = Extract3DNeighborhoods(Sgrain, nbOffsets, Recon.NUM_CORES);
                
                % Convert the grain list to a uint16. It takes up much less space
                GrainList = uint16(GrainList);

                % Get the number of grains
                NumGrains = sum(GrainList > 0, 2);

                % Truncate the grain list, this gets rid of zeros.
                GrainList = GrainList(:, 1:max(NumGrains));

                [Xnnidx] = flann_search(Recon.Exemplar_Index{ExIndex}, [NB_queries NumGrains*1000]', 1, Recon.Exemplar_Params{ExIndex});
            else
                [Xnnidx] = flann_search(Recon.Exemplar_Index{ExIndex}, NB_queries', 1, Recon.Exemplar_Params{ExIndex});
            end
        end

        % Now, mark any neighborhoods that are on the edges as a -1, this
        % will tell the code later to ignore these guys in our calculations
        Xnnidx(isPeriodic) = -1;
        
        % We need to reshape the search results into the 4D array that the
        % rest of our code expects. We need to do the permute because of
        % how we arranged the neighborhoods when extracting them from the
        % 3D image. We did it in the opposite of the MATLAB convention for
        % some reason, maybe because they are extracted in C++ code.
        Recon.NNB_Table(:,:,:,ExIndex) = permute(reshape(Xnnidx, size(S)), [3 2 1]);
        
    end

end
