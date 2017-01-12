function [Recon] = DestroyRecon(Recon)
%DestroyRecon - This function deallocates any memory objects 
% created when calling SetupRecon. Every call to SetupRecon must have a
% corresponding call to DestroyRecon. If not, your program will have memory
% leaks! Do not try to use a Recon struct after it has been destroyed,
% unexpected behaviours will occur!
%
% Syntax:  DestroyRecon(Recon)
%
% Inputs:
%    Recon - The struct that was created by SetupRecon.
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%

% Author: David Turner
% email: davidt0x@gmail.com
% December 2016

%------------- BEGIN CODE --------------

% Clears up the memory associated with any ANN index
if(isfield(Recon, 'Exemplar_Index'))
    
    % If we are working with FLANN indices
    if(strcmp(Recon.ANN_ALGO, 'FLANN'))
        
        % Loop over each of the exemplar images indices.
        for ii=1:size(Recon.Exemplar_Index, 1)
            
            % If the index has not been freed already.
            if(Recon.Exemplar_Index{ii} ~= 0)
                fprintf(1, 'Freeing FLANN Index\n');
                flann_free_index(Recon.Exemplar_Index{ii}); 
            end
            
            % Mark the index as free by setting the pointer equal to 0.
            Recon.Exemplar_Index{ii} = 0;
        end
    
    end
    
    % If we are working with PatchTable indices
    if(strcmp(Recon.ANN_ALGO, 'PatchTable'))
        
        % Loop over each of the exemplar images indices.
        for ii=1:size(Recon.Exemplar_Index, 1)
            
            % If the index has not been freed already.
            if(Recon.Exemplar_Index{ii} ~= 0)
                fprintf(1, 'Freeing PatchTable Index\n');
%                flann_free_index(Recon.Exemplar_Index{ii}); 
            end
            
            % Mark the index as free by setting the pointer equal to 0.
            Recon.Exemplar_Index{ii} = 0;
        end
    
    end
    
end