function [] = DestroyReconMultiRes(ReconHierarchy)
%DestroyReconMultiRes - This function deallocates any memory objects 
% created when calling SetupReconMultiRes. Every call to SetupReconMultiRes must have a
% corresponding call to DestroyReconMultRes. If not, your program will have memory
% leaks! Do not try to use a Recon struct after it has been destroyed,
% unexpected behaviours will occur!
%
% Syntax:  DestroyReconMultiRes(ReconHierarchy)
%
% Inputs:
%    ReconHierarchy - The struct that was created by SetupReconMultiRes
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%

% Author: David Turner
% email: davidt0x@gmail.com
% December 2016

%------------- BEGIN CODE --------------

    if(isfield(ReconHierarchy, 'ReconObjects') && iscell(ReconHierarchy.ReconObjects) && ~isempty(ReconHierarchy.ReconObjects)) 
        % Loop trough the reconstruction objects and call destory on each one.
        for ii=1:length(ReconHierarchy.ReconObjects)
            DestroyRecon(ReconHierarchy.ReconObjects{ii});
        end
    end
    
    ReconHierarchy.ReconObjects = cell.empty(0);
    
end
