% This function is called at each iteration to SolidOptimization
% to plot the current results. This is called to plot a grain map
function [] = PlotIterationPolycrystal(Recon, S, Snew, useWeights, H)

    PAD = round(Recon.NB_SIZE/2);
    
    Snew = Snew(PAD:(size(Snew,1)-PAD), PAD:(size(Snew,2)-PAD), PAD:(size(Snew,3)-PAD));

    rng(1);
    cmap = rand(3000,3);

    subplot(1,2,1)
    MID = ceil(size(Snew)/2);
    slice(Snew, [1 size(Snew,1)], [1 size(Snew,2)], [1 size(Snew,3)]); 
    axis equal; axis tight; shading flat; box on; 
    %colormap('gray'); caxis([0 1]); 
    colormap(cmap); caxis([1 3000]);
    freezeColors;
    axis off;
    colorbar('off')
    subplot(1,2,2);
    slice(Snew, MID(1), MID(2), MID(3)); 
    axis equal; axis tight; shading flat; box on; 
    colorbar('off')
    %colormap('gray'); caxis([0 1]);
    colormap(cmap); caxis([1 3000]);
    freezeColors;
    axis off;

%    subplot(2,2,4);
%    T = Recon.SourceTable(:, :, :, 1);  
%    T = T(:)+1;
%    P = zeros(size(Recon.E{1}));
%    for zz=1:length(T)
%        P(T(zz)) = P(T(zz)) + 1;
%    end
%    P(P == 0) = NaN;
    
%    if(useWeights)
%        H = H(1:end-Recon.NB_SIZE, 1:end-Recon.NB_SIZE);
%        H(H == 0) = NaN;
%        pcolor(flipud(H)); shading flat; colorbar; colormap('jet');
%    else
%        imagesc(P); colormap('jet');
%    end


   drawnow;

end
