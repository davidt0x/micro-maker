function [isValid, errMsg] = CheckReconParams(recon_params)
%CheckReconParams - This function checks whether the reconstruction
%parameters are valid. 
%
% Syntax:  [isValid errMsg] = CheckReconParams(recon_params)
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
  
    % Now lets see if they have the neighborhood sizes correct
    if(~isfield(recon_params, 'NB_SIZE') || length(recon_params.NB_SIZE) ~= 1)
        errMsg = 'NB_SIZE must exist and it must be a scalar. All neighborhoods are square.';
        return;  
    end

    % All neighborhood sizes must be an odd numebr of pixels and be
    % at least 3 pixels.
    if(~mod(recon_params.NB_SIZE,2) || recon_params.NB_SIZE < 3)
        errMsg = 'NB_SIZE must be odd and greater that 1';
        return;
    end
        
    if(~isfield(recon_params, 'RECON_SIZE') || length(recon_params.RECON_SIZE) ~= 3 || ... 
            any(recon_params.RECON_SIZE) < 1)
        errMsg = 'RECON_SIZE must exist and its length must be 3 and all values must be greater than 0';
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
       
        % Each element of the cell array must be a double matrix
        % with values ranging between 0 and 1 
        if(~ismatrix(E))
            errMsg = 'Elements of EXEMPLARS must be 2D matrices';
            return
        end
        
        % If we have a cell array of images then check each one
        [isOK, exErrMsg] = IsValidExemplar(E);
        if(~isOK)
            errMsg = exErrMsg;
            return
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