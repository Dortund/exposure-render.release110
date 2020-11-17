function stat = CompareImagesWithName(gt, img, maxvalue, name)

if ( nargin < 3 )
    maxvalue = 255.0;
end

diffsquared = (gt - img).^2;

%imwrite(((gt - img) * 255 + 255) ./ 2, name);
imwrite(uint8(diffsquared*255*100), name);

% sum of squared differences
stat.ssd     = sum(  sum( diffsquared ) );
stat.ssdsum  = sum(  stat.ssd );
% display( sprintf('SSD for first  channel  is %15.9f', stat.ssd(1) ) );
% if (length(stat.ssd) > 1)
% display( sprintf('SSD for second channel  is %15.9f', stat.ssd(2) ) );
% end
% if (length(stat.ssd) > 2)
% display( sprintf('SSD for third  channel  is %15.9f', stat.ssd(3) ) );
% end
% display( sprintf('SSD for all    channels is %15.9f\n', stat.ssdsum ) );


% mean squared error
stat.mse     = mean( mean( diffsquared ) );
stat.msesum  = sum(  stat.mse );
% display( sprintf('MSE for first  channel  is %15.9f', stat.mse(1) ) );
% if (length(stat.mse) > 1)
% display( sprintf('MSE for second channel  is %15.9f', stat.mse(2) ) );
% end
% if (length(stat.mse) > 2)
% display( sprintf('MSE for third  channel  is %15.9f', stat.mse(3) ) );
% end
% display( sprintf('MSE for all    channels is %15.9f\n', stat.msesum ) );


% root mean squared error
stat.rmse     = sqrt( stat.mse );
stat.rmsesum  = sqrt( stat.msesum );
% display( sprintf('RMSE for first  channel  is %15.9f', stat.rmse(1) ) );
% if (length(stat.rmse) > 1)
% display( sprintf('RMSE for second channel  is %15.9f', stat.rmse(2) ) );
% end
% if (length(stat.rmse) > 2)
% display( sprintf('RMSE for third  channel  is %15.9f', stat.rmse(3) ) );
% end
%display( sprintf('RMSE for all    channels is %15.9f\n', stat.rmsesum ) );


% peak signal-to-noise-ratio
stat.psnr     = 20 * log10( maxvalue / stat.rmse );
stat.psnrsum  = 20 * log10( maxvalue / stat.rmsesum );
% display( sprintf('PSNR for first  channel  is %15.9f', stat.psnr(1) ) );
% if (length(stat.psnr) > 1)
% display( sprintf('PSNR for second channel  is %15.9f', stat.psnr(2) ) );
% end
% if (length(stat.psnr) > 2)
% display( sprintf('PSNR for third  channel  is %15.9f', stat.psnr(3) ) );
% end
% display( sprintf('PSNR for all    channels is %15.9f\n', stat.psnrsum ) );


% Structural SIMilarity (SSIM)
[stat.ssim(1),stat.ssim_index(1,:,:)] = ssim( gt(:,:,1), img(:,:,1));%, [0.01 0.03], fspecial('gaussian', 11, 1.5), 1.0);
if (length(stat.ssim) > 1)
[stat.ssim(2),stat.ssim_index(2,:,:)] = ssim( gt(:,:,2), img(:,:,2));%, [0.01 0.03], fspecial('gaussian', 11, 1.5), 1.0);
end
if (length(stat.ssim) > 2)
[stat.ssim(3),stat.ssim_index(3,:,:)] = ssim( gt(:,:,3), img(:,:,3));%, [0.01 0.03], fspecial('gaussian', 11, 1.5), 1.0);
end
stat.ssimsum  = sum(  stat.ssim );
stat.ssimmean = mean( stat.ssim );
% display( sprintf('SSIM        for first  channel  is %15.9f', stat.ssim(1) ) );
% if (length(stat.ssim) > 1)
% display( sprintf('SSIM        for second channel  is %15.9f', stat.ssim(2) ) );
% end
% if (length(stat.ssim) > 2)
% display( sprintf('SSIM        for third  channel  is %15.9f', stat.ssim(3) ) );
% end
% display( sprintf('Sum of SSIM for all    channels is %15.9f', stat.ssimsum ) );
% display( sprintf('Mean SSIM   for all    channels is %15.9f\n', stat.ssimmean ) );

% end function


