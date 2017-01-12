% This function is called at each iteration to SolidOptimization
% to plot the current results.
function [] = PlotIteration(Recon, S, Snew, useWeights)

    % Check to make sure that Snew is not all of one value. This
    % will cause slice to error
    if(max(Snew(:)) ~= min(Snew(:)))

        % In the first subplot, show the outer faces of the cubic volume in
        % greyscale
        subplot(2,2,1);
        MID = ceil(size(Snew)/2);
        slice(Snew, [1 size(Snew,1)], [1 size(Snew,2)], [1 size(Snew,3)]); 
        axis equal; axis tight; shading flat; box on; 
        colormap('gray'); caxis([0 1]); 
        freezeColors;
        axis off;
        
        % In the second subplot, show the three intersecting orthgonal 
        % medial slices of the cubic volume in grey scale. 
        subplot(2,2,2);
        slice(Snew, MID(1), MID(2), MID(3)); 
        axis equal; axis tight; shading flat; box on; 
        colormap('gray'); caxis([0 1]);
        freezeColors;
        axis off;

    end

    % In the third subplot, show the histogram of values for the
    % reconstruction and one of the exemplars for comparison.
    subplot(2,2,3);
    Shist = histc(Snew(:), 0:(1/16):1);
    Ehist = histc(Recon.EXEMPLARS{1}(:), 0:(1/16):1);
    Shist = Shist ./ sum(Shist);
    Ehist = Ehist ./ sum(Ehist);
   
    hold off;
    plot(0:(1/16):1, Shist, '-r', 'DisplayName', 'Sample');
    hold on;
    plot(0:(1/16):1, Ehist, '-b', 'DisplayName', 'Exemplar');
    title('Value Histogram');
    
    subplot(2,2,4);
    T = Recon.SourceTable(:, :, :, 1);  
    T = T(:)+1;
    P = zeros(size(Recon.EXEMPLARS{1}));
    for zz=1:length(T)
        P(T(zz)) = P(T(zz)) + 1;
    end

    P(P == 0) = NaN;
    %imagesc(weight_table); colorbar; colormap('jet');
    
    if(useWeights)
        H = Recon.NBHoodHist{1};
        H = H(1:end-Recon.NB_SIZE, 1:end-Recon.NB_SIZE);
        H(H == 0) = NaN;
        pcolor(flipud(H)); shading flat; colorbar; colormap('jet');
        %H2 = H(:);
        %H2 = H2(~isnan(H2));
        %caxis([mean(H2(:))-3*std(H2(:)) mean(H2(:))+3*std(H2(:))])
        title('Exemplar #1 Neighborhood Histogram');
    else
        pcolor(flipud(P)); shading flat; colormap('jet'); colorbar; axis equal; axis tight;
        title('Exemplar #1 Pixel Histogram');
    end

    drawnow;

end
