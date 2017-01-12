% Lets do a quick check to see if the mex files have been built. If not, we
% will die with and error.
d = dir('../cpp/bin');

if(sum(~cellfun(@isempty,strfind({d.name}, 'mex'))) == 0)
   error('MEX code does not appear to be built. Try running build_mex first'); 
else
    addpath('../cpp/bin');
end

% Check if we can find FLANN on the path, if not, error out
if(exist('flann_build_index') == 0)
    error('FLANN does not appear on the path. Modify path to add it.')
end

%addpath('../cpp/PatchTableMEX');