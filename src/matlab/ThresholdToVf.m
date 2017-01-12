function [M2] = ThresholdToVf(M, vf_to_match);

    vf_match = @(M, vf) mean(double(M(:) > vf));

    lvls = 0:0.001:1;
    vfs = zeros(1, length(lvls));
    ii = 1;
    for vf=lvls
        vfs(ii) = vf_match(M, vf);
        ii = ii + 1;
    end

    [minV minI] = min(abs(vfs - vf_to_match));

    M2 = double(M > lvls(minI));

end
