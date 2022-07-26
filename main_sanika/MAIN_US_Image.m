clc
clear all
close all
warning off

Fn = dir('.\UltraSound Images\*.*');
inc = 1;
for I = 3 : 7%length(Fn)
    
    if ~strcmp(Fn(I).name,'.') == 1 && ~strcmp(Fn(I).name,'..') == 1 && ~strcmp(Fn(I).name,'thumbs.db') == 1
        
        %% =============== Image Reading ==============================
        
        Img = imresize(imread(['.\UltraSound Images\' Fn(I).name ]),[256 256]);
        
        % convert Rgb to gray
        
        if size(Img,3) > 1
            gImg = rgb2gray(Img);
        else
            gImg  = Img;
        end
        
        figure;imshow(gImg);title('Input Image')
        
        
        
        % Apply Speckle Noise
        
        NImg1 =  imnoise(gImg,'Speckle', 0.1);
        [Noisy_psnr Noisy_snr] = psnr(im2double(NImg1),im2double(gImg));
        NoisyPsnrValue(I,1) = Noisy_psnr;
        NoisySnrValue(I,1) = Noisy_snr;
        figure;imshow(NImg1);title('Noisy Image')
        
        %% Apply Median filter
        
        
        F_img2 = medfilt2(NImg1);
        
        % Apply modified wiener filter
        
        F_img1 = MMWF_2D(im2double(F_img2),5);
        
        % Apply Bilateral Filter
        
        F_img = im2double(F_img1);
        w     = 2;       % bilateral filter half-width
        sigma = [3 0.1]; % bilateral filter standard deviations
        F_img(F_img<0) = 0; F_img(F_img>1) = 1;
        
        BImg = bfilter2(double(F_img),w,sigma);
        
        
        %% ============= Bivariate Shrinkage =====================
        
        bivrnt_DeImg = bishrink(im2double(BImg),im2double(gImg),0.0009);
        
        %% ============= ETV denoising process =====================
        
        
        lambda = 0.3; mu = 0.4; nu = 0.5; iter = 40; %wdc good psnr with funSSTV
        ETV_DeImg = funETV(im2double(BImg),iter,lambda,mu,nu);
        
        
        
        %% ============= Image Fusion using DTCWT ===================
        
        Fsd_Img = DTCWT_IDTCWT(double(bivrnt_DeImg),double(ETV_DeImg));
        
        %% ============= Debluring Using Optimized Richard Luzy Technique =======
        
        % WithOut Optimization
        
        PSF = fspecial('gaussian',2,7);
        RL_Output = deconvlucy(Fsd_Img,PSF);
        figure;imshow(RL_Output,[]);title('RL Denoised Image')
        [RL_psnr, RL_snr] = psnr(im2double(RL_Output), im2double(gImg));
        RL_Psnr_Value(I,1) = RL_psnr;
        RL_Snr_Value(I,1) = RL_snr;
        
        % With Optimization
        
        [Fitness RL_Op_Output PSF] = KH_Optimization_RL1(im2double(Fsd_Img),im2double(gImg));
        figure;imshow(RL_Op_Output,[]);title('KH RL Opt Denoised Image')
        [RL_opt_psnr, RL_opt_snr] = psnr(im2double(RL_Op_Output), im2double(gImg));
        OP_RL_Psnr_Value(I,1) = RL_opt_psnr;
        OP_RL_Snr_Value(I,1) = RL_opt_snr;
        
        Res = [Noisy_psnr Noisy_snr RL_psnr, RL_snr RL_opt_psnr, RL_opt_snr];
        Result(inc,:) = Res;
        inc = inc + 1;
    end
end