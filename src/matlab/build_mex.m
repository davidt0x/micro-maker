% This script builds the C++ mex portion of this software. If mex -setup
% has been configured correctly then it should build these files on most
% systems without issue, as they do not depend on external libraries. I
% have currently tested this build process under Windows 10 with Visual
% Studio 2013 as the compiler. Minor changes might need to be made to
% compile under GNU. 

% Clear any loaded mex files so we don't have problems overwriting the mex
% files
clear mex;

fprintf(1, 'Building OptimizeSolidMEX\n');
mex -outdir ../cpp/bin/ ../cpp/OptimizeSolidMEX.cpp

fprintf(1, 'Building Extract3DNeighborhoods\n');
mex -outdir ../cpp/bin/ ../cpp/Extract3DNeighborhoods.cpp

fprintf(1, 'Building GetAllNHoodsMEX\n');
mex -outdir ../cpp/bin/ ../cpp/GetAllNHoodsMEX.cpp
