function [S, Recon, S_iters] = SolidOptimization(Recon, NUM_ITERATIONS, start, useWeights)

    % Set the default arguments

    % If the user didn't specify whether we should use neighborhood
    % weighing, default to no.
    if(nargin < 4)
        useWeights = 0;
    end

    % Lets set the iteration variable to the starting guess
    S = start;

    % If we are using neighborhoods weights, we might need to initialize 
    % them first. If they exit already, don't worry about it.
    if(useWeights)
        if(~isfield(Recon, 'NBWeights'))
            fprintf(1, 'Initializing Neighborhood Weights\n');
            for ii=1:length(Recon.EXEMPLARS)
                Recon.NBWeights{ii} = ones(1, size(Recon.Exemplar_NBs{ii},1)) * Recon.NB_SIZE;
            end
        end
    end
    
    % If the user has asked to return a copy of each iteration's structure
    % then we need to allocate it.
    if(nargout > 2)
        S_iters = zeros([size(S) NUM_ITERATIONS+1]);
        S_iters(:, :, :, 1) = S;
    end
    
    elapsed_search = 0;
    elapsed_opt = 0;

    for ii=1:NUM_ITERATIONS
        fprintf(1, 'Iteration: %d\n', ii);

        tic;
        fprintf(1, '    Finding nearest neighbors ... ');
        Recon = SearchNNB(S, Recon, useWeights);
        elapsed_search = elapsed_search + toc;
        fprintf(1, 'Done: %f seconds\n', toc);

        tic;
        fprintf(1, '    Optimizing solid ... ');
        [Snew, Recon.TexelWeights, Recon.SourceTable] = OptimizeSolidMEX(S, Recon);

        % For each exemplar, lets calculate a histogram that keeps track
        % of how often that neighborhood has been matched. This will give
        % us and idea of how much of the exemplar is being used in each
        % iteration
        Recon.NBHoodHist = cell(length(Recon.EXEMPLARS),1);
        for pp=1:size(Recon.EXEMPLARS, 1)
            % Calculate a histogram that counts many times each
            % neighborhood has been selected as the best matching
            % neighborhood.
            H = zeros(size(Recon.EXEMPLARS{pp})); 
            NBs = Recon.NNB_Table(:, :, :, pp); 
            NBs = NBs(:);
            NBs = NBs(NBs > 0);
            for kk=1:length(NBs) 
                xy = Recon.NB_ExemplarLookup{pp}(NBs(kk), :); 
                H(xy(1), xy(2)) = H(xy(1), xy(2)) + 1; 
            end
            H = H ./ sum(H(:));
            
            Recon.NBHoodHist{pp} = H;
        end
        
        % Update the neighborhood weights
        if(useWeights)
            for pp=1:size(Recon.EXEMPLARS, 1)
                
                % If every neighborhood is used equally, then we would
                % expect to find this value for each neighborhood 
                expected_weight = (length(NBs) / size(Recon.Exemplar_NBs{pp}, 1)) / length(NBs);
                
                % Calculate how far off we are for each neighborhood
                diff_weights = Recon.NBHoodHist{pp} - expected_weight;
                  
                % This scale factor is a fudge factor. Basically, it
                % controls how much the weights will change. The larger the
                % value, the quicker neighborhoods will get downweighted
                % each iteration. Remember, the neighborhood weights will
                % keep increasing for oversampled neighborhoods each
                % iteration if the image wont change. This allows us to get
                % out of local minima if we wait enough iterations.
                ScaleFactor = 10*Recon.NB_SIZE*Recon.NB_SIZE;

                % Accumulate the changes to the weight table 
                weight_table = zeros(size(diff_weights));
                for kk=1:length(Recon.NBWeights{pp})
                    xy = Recon.NB_ExemplarLookup{pp}(kk, :);
                    Recon.NBWeights{pp}(kk) = Recon.NBWeights{pp}(kk) + diff_weights(xy(1), xy(2)) * ScaleFactor;
                    if(Recon.NBWeights{pp}(kk) < 0)
                        Recon.NBWeights{pp}(kk) = 0;
                    end
                    weight_table(xy(1), xy(2)) = Recon.NBWeights{pp}(kk);
                end
            end
        end

        % If we are running a polycrsytal reconstruction
        if(isfield(Recon, 'ANG_MAP'))
            PlotIterationPolycrystal(Recon, S, Snew, useWeights, ANG_MAP);            
        else
            PlotIteration(Recon, S, Snew, useWeights);
        end
       
        % Store a copy of the current structure every iteration if the user
        % wants it.
        if(nargout > 2)
            S_iters(:, :, :, ii+1) = Snew;
        end

        % Calculate the percentage of the image that has changed
        perVoxelPercentChange = abs(S(:) - Snew(:)) ./ S(:); 
        
        % If we have a NaN, that means 0 / 0, which we will say is no
        % relative change.
        perVoxelPercentChange(isnan(perVoxelPercentChange)) = 0;
        
        % If we have Inf, then we have something like x/0. Which we will
        % just say is 1.
        perVoxelPercentChange(isinf(perVoxelPercentChange)) = 1;
        
        % Calculate the total percentage change over the entire image.
        % This is the relative percent change divided by the total possible
        % change.
        percentChange = 100*(sum(perVoxelPercentChange) / numel(S));
          
        % It converged!
        if( ii == NUM_ITERATIONS )
            fprintf(1, 'Done.\nIt converged!\n');
            break;
        else
            S = Snew;
        end

        fprintf(1, 'Done: %f, vf=%f seconds, percentChange=%f\n', toc, mean(S(:)), percentChange);
        elapsed_opt = elapsed_opt + toc;

    end
    fprintf(1, '   Iterations: %d   SearchTime: %f secs   OptimTime: %f secs\n', ii, elapsed_search, elapsed_opt);

    % Resize the iterations matrix if they want it. We might not need it 
    % all if the reconstruciton converged earlier that NUM_ITERATIONS
    if(nargout > 2)
        S_iters = S_iters(:, :, :, 1:ii);
    end

end
