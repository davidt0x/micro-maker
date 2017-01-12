function [S_star, Recon] = Reconstruct(Recon, max_iterations, startm)
% Reconstruct - This function runs the reconstruction algorithm using the
% reconstruction struct passed to it. The Recon object must be created with
% either SetupRecon or SetupReconMultiRes.
%
% Syntax:  [S_star, Recon] = Reconstruct(Recon, max_iterations, startm)
%
% Inputs:
%    Recon - A reconstruction struct created by either SetupRecon or
%    SetupReconMultiRes.
%
%    max_iterations - Either a vector of length Recon.NUM_LEVELS or a
%    scalar specifiying the number of iterations to run the reconstruction
%    at each level. The levels are ordered from highest resolution to
%    lowest resolution. It is wise to choose a low number of resolutions
%    for higher resolutions because they take much longer. For example, for
%    a three level reconstruction something like [5, 50, 100] might work
%    well. This is very problem dependent thought as some reconstructions
%    take longer to converge than others.
%
%    startm - The starting guess microstructure to work with. This must be
%    the size of the lowest level resolution in your reconstruction
%    hierarchy. That is, Recon.ReconObjects{Recon.NUM_LEVELS}.RECON_SIZE.
%
% Outputs:
%    RS - The reconstruction struct containing all data structures
%    needed for performing a reconstruction optimization.
%
% Example: 
%    RS = SetupRecon(params);
%
% Other m-files required: 
% MAT-files required: ThresholdToVf.m PlotIteration.m. PlotIterationGrainMap.m
%
% See also: SetupRecon SetupReconMultiRes 

% Author: David Turner
% Email: davidt0x@gmail.com 
% Website: https://github.com/davidt0x
% December 2016

%------------- BEGIN CODE --------------

    % We need to check whether this is hierarchical reconstruction or just
    % a single level reconstruction.
    if(isfield(Recon, 'NUM_LEVELS'))
        FULL_RECON_SIZE = Recon.ReconObjects{1}.RECON_SIZE;
        START_RECON_SIZE = Recon.ReconObjects{Recon.NUM_LEVELS}.RECON_SIZE;
    else
        FULL_RECON_SIZE = Recon.RECON_SIZE;
        START_RECON_SIZE = FULL_RECON_SIZE;
    end

    % Setup all the default arguments
    
    if(nargin < 2)
        if(isfield(Recon, 'NUM_LEVELS'))
            max_iterations = ones(Recon.NUM_LEVELS, 1) * 100;
            max_iterations(1) = 10;
        else
            max_iterations = 100;
        end
    end

    
    % If the user has not passed in a starting guess then create one
    % from threhsolded random noise.
    if(nargin < 3)
       startm = double(rand(START_RECON_SIZE) > 0.5);
    else
        % Make sure the dimensions are correct
        if( any(size(startm) ~= START_RECON_SIZE) )
            error('Starting guess microstructure has incorrect dimensions. Must be [%d %d %d]', ... 
                START_RECON_SIZE(1), START_RECON_SIZE(2), START_RECON_SIZE(3));
        end
    end
    
    % If the user did not pass in whether to use neighborhood search
    % waiting then assume the default (no)
    if(~isfield(Recon, 'UseNeighborhoodSearchWeights'))
        useWeights = false;
    else
        useWeights = Recon.UseNeighborhoodSearchWeights;
    end

    % Lets initialize the loop reconstruciton variable witht the starting
    % guess.
    S_star = startm;
       
    % We need to check whether this is hierarchical reconstruction or just
    % a single level reconstruction.
    if(isfield(Recon, 'NUM_LEVELS'))
        
        % If we are using neighborhood weighting, then we need to remove
        % any weights from previous runs. This will trigger the first call
        % to SolidOptimization to re-initialize them.
        if(useWeights)
            for ii=1:Recon.NUM_LEVELS
                if(isfield(Recon.ReconObjects{ii}, 'NBWeights'))
                    Recon.ReconObjects{ii} = rmfield(Recon.ReconObjects{ii}, 'NBWeights');
                end
            end
        end

        
        % Lets initialize the loop variables
        ll = 1;
        
        % We will be plotting a lot of figures during the optimization
        % process, close anything that is open.
        close all;

        % Lets start the reconstruction, we loop through the reconstruction
        % objects from last to first because this is from lowest resolution
        % to highest.
        for ii=Recon.NUM_LEVELS:-1:1

            % Don't use neighborhood weights at the highest resolution, it takes too long.
            if(ii == 1)
                useWeights = 0;
            end

            % Make sure to assign the return reconstruction object because
            % this functions is expected to produce side effects on the
            % object.
            [S_star Recon.ReconObjects{ii}] = SolidOptimization(Recon.ReconObjects{ii}, ...
                                                max_iterations(ii), S_star, useWeights);
   
            
            % If we are not at the final resolution, we need to do some
            % things before going into the next resolution. Basically, we
            % need upsample three things, the current reconstruction
            % results, the weight table for texels, and the weight table
            % for neighborhoods. Otherwise, just go to the next iteration
            % (break)
            if(ii == 1)
                break;
            end

            % Upsample the volume. Get the next iterations size. 
            nxSize = Recon.ReconObjects{ii-1}.RECON_SIZE;

            % Upsample the current results, we will use simple nearest
            % neighbor interpolation to avoid blurring. Hopefully, it
            % produces good results.
            S_star = resize(S_star, nxSize, 'nearest');

            % If we are going to the final level next then we need to threshold
            % the input. Lets choose a threshold such that the volume fraction
            % matches the exemplar.
            if(ii == 2)

                % Compute the average volume fraction accross all
                % exemplars.
                vf_to_match = 0;
                for zz=1:length(Recon.EXEMPLARS)
                    T = Recon.EXEMPLARS{zz};
                    vf_to_match = vf_to_match + mean(T(:));
                end
                vf_to_match = vf_to_match / length(Recon.EXEMPLARS);

                % Threshold the image so the volume fraction is as
                % close as we can get.
                S_star = ThresholdToVf(S_star, vf_to_match);

            end

            % Resize the texel weights so we can use them to initialize
            % the next level.
            for pp=1:length(Recon.ReconObjects{ii}.EXEMPLARS)
                Recon.ReconObjects{ii-1}.TexelWeights{pp} = ...
                    resize(Recon.ReconObjects{ii}.TexelWeights{pp}, ... 
                           size(Recon.ReconObjects{ii-1}.EXEMPLARS{pp}), 'nearest');
            end

            % Ok, if we have a neighborhood weighting table, resize it 
            % as well and use it as the input into the next level. Only 
            % do this if we are not going into the final resolution 
            % because we will not use neighborhood weighting at the 
            % final resolution.
            if(useWeights && ii > 2)

                % The neighborhood weight table is not stored like an image but instead
                % a list. We need to make it an image so we can resize\interpolate spatially
                for pp=1:length(Recon.ReconObjects{ii}.EXEMPLARS)
                    weight_table = zeros( size(Recon.ReconObjects{ii}.EXEMPLARS{pp})-Recon.ReconObjects{ii}.NB_SIZE+1 );

                    % Reformat NBWeights{pp} into an image.
                    for kk=1:length(Recon.ReconObjects{ii}.NBWeights{pp})
                        % Get the x,y pixel location of this neighborhood in the exemplar
                        xy = Recon.ReconObjects{ii}.NB_ExemplarLookup{pp}(kk, :);

                        weight_table(xy(1), xy(2)) = Recon.ReconObjects{ii}.NBWeights{pp}(kk);
                    end

                    % Upsample the weight_table with simple nearest neighbor interpolation
                    nx_Ex_Size = size(Recon.ReconObjects{ii-1}.EXEMPLARS{pp})-Recon.ReconObjects{ii}.NB_SIZE+1;
                    weight_table = resize(weight_table, nx_Ex_Size, 'nearest');

                    % Reformat it as a list.
                    Recon.ReconObjects{ii-1}.NBWeights{pp} = zeros(1, size(Recon.ReconObjects{ii-1}.Exemplar_NBs{pp},1));
                    for kk=1:length(Recon.ReconObjects{ii-1}.NBWeights{pp})
                        xy = Recon.ReconObjects{ii-1}.NB_ExemplarLookup{pp}(kk, :);
                        Recon.ReconObjects{ii-1}.NBWeights{pp}(1, kk) = weight_table(xy(1), xy(2));   
                    end

                end

            end
            
            ll = ll + 1;
        
        end % For Each Level Of Reconstruction

    else
        % We have just a single level. We can just call SolidOptimization on
        % it.
        [S_star, Recon] = SolidOptimization(Recon, max_iterations, S_star, useWeights);
    end

end