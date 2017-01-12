function [isValid, errMsg] = CheckReconParamsMultiRes(recon_params)
%CheckReconParamsMultiRes - This function checks whether the reconstruction
%parameters are valid for a multi-resolution reconstruction. 
%
% Syntax:  [isValid errMsg] = CheckReconParamsMultiRes(recon_params)
%
% Inputs:
%    recon_params - The reconstruction parameters specified as a struct
%
% Outputs:
%    isValid - True if everything looks valid, false otherwise.
%    errMsg - If it isn't valid, this is the reason.
%
% Other m-files required: none
% Subfunctions: IsValidExemplar(Image)
% MAT-files required: none
%

% Author: David Turner
% email: davidt0x@gmail.com
% December 2016

%------------- BEGIN CODE --------------

    isValid = false;
    
    % Lets check to see if the number of levels is specified correctly
    if(~isfield(recon_params, 'NUM_LEVELS') || recon_params.NUM_LEVELS < 1)
        errMsg = 'NUM_LEVELS must exist and be greater than 0';
        return;  
    end

    % Now lets see if they have the neighborhood sizes correct
    if(~isfield(recon_params, 'NB_SIZES') || length(recon_params.NB_SIZES) ~= recon_params.NUM_LEVELS)
        errMsg = 'NB_SIZES must exist and its lenght must be equal to NUM_LEVELS';
        return;  
    end

    % All neighborhood sizes must be an odd numebr of pixels and be
    % at least 3 pixels.
    if(any(~mod(recon_params.NB_SIZES,2)) || any(recon_params.NB_SIZES < 3))
        errMsg = 'NB_SIZES must all be odd and greater that 1';
        return;
    end
        
    if(~isfield(recon_params, 'FULL_RECON_SIZE') || length(recon_params.FULL_RECON_SIZE) ~= 3 || ... 
            any(recon_params.FULL_RECON_SIZE) < 1)
        errMsg = 'FULL_RECON_SIZE must exist and its lenght must be 3 and all values must be greater than 0';
        return;  
    end

    if(~isfield(recon_params, 'EXEMPLARS') || ~iscell(recon_params.EXEMPLARS) || isempty(recon_params.EXEMPLARS))
        errMsg = 'EXAMPLARS must be a cell array containing one image per element';
        return;  
    end

    function [isOK, errMsg] = IsValidExemplar(Image)
        isOK = false;
        errMsg = '';
        if(~isa(Image, 'double'))
           errMsg = 'Each exemplar image must be a matrix of type double';
           return;
        end      
        if(max(Image(:)) > 1 || min(Image(:)) < 0)
           errMsg = 'Exemplar image pixels must all range between 0 and 1';
           return;
        end    
        isOK = true;
    end
    
    % Lets check that the exemplar images are the appropriate form.
    for ii=1:length(recon_params.EXEMPLARS)
        E = recon_params.EXEMPLARS{ii};
       
        % Each element of the cell array must be either a double matrix
        % with values ranging between 0 and 1 or another cell array that
        % contains images like these.
        if(~ismatrix(E) && ~iscell(E))
            errMsg = 'Elements of EXEMPLARS must be either 2D matrices or cell arrays';
            return
        end
        
        % If we have a cell array of images then check each one
        if(iscell(E))
          for jj=1:length(E)
            [isOK exErrMsg] = IsValidExemplar(E{jj});
            if(~isOK)
                errMsg = exErrMsg;
                return
            end 
          end
        % If this is just an image then check it.
        else
            [isOK exErrMsg] = IsValidExemplar(E);
            if(~isOK)
                errMsg = exErrMsg;
                return
            end    
        end
    end
    
    if(~isfield(recon_params, 'NB_INDICES') || ... 
            length(recon_params.NB_INDICES) ~= length(recon_params.EXEMPLARS) || ...
            any(recon_params.NB_INDICES < 1 | recon_params.NB_INDICES > 9))
       errMsg = 'NB_INDICES must be a vector of the same length as EXEMPLARS and each element must be between 1 and 9'; 
       return;
    end
    
    if(isfield(recon_params, 'ANN_ALGO') && ... 
       ~strcmp(recon_params.ANN_ALGO, 'FLANN') && ... 
       ~strcmp(recon_params.ANN_ALGO, 'PatchTable'))
        errMsg = 'ANN_ALGO must be set to either PatchTable or FLANN';
        return;
    end
  
    errMsg = '';
    isValid = true;
end