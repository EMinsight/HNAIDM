function [ lithoImage, Intensity ] = Calculate2DAerialImage(source, mask, projector, lithoImage, receipe, numerics)

%% 数据初始化
maskNf = numerics.SampleNumber_Mask_X;  % 默认81
maskNg = numerics.SampleNumber_Mask_Y;  % 默认81
waferNf = numerics.SampleNumber_Wafer_X; % 默认81
waferNg = numerics.SampleNumber_Wafer_Y; % 默认81

source.PntNum = numerics.SampleNumber_Source;
[ sourceData, ~ ] =  source.Calc_SourceSimple();
weight = sum( sourceData.Value ); 
wavelength = source.Wavelength; % 193
NA = projector.NA;
M = projector.Reduction;
if strcmpi(projector.LensType,'Immersion')
    indexImage = projector.Index_ImmersionLiquid;
elseif strcmpi(projector.LensType,'Dry')
    indexImage = 1;   % 1.44
    if NA>=1
        error('Wrong NA!');
    end
else
    error('Unsupported Lens Type')
end

% 掩模频谱
mask.Nf = maskNf;
mask.Ng = maskNg;
[ spectrum, mask_fs, mask_gs ] = mask.CalculateMaskSpectrum(projector ,source);

SimulationRange = receipe.FocusRange - receipe.Focus;

%% 开始计算
Intensity = zeros(waferNf,waferNg,length(SimulationRange));
indexK = indexImage;

Orientation = mask.Orientation;

for iFocus = 1:length(SimulationRange)

    [ mask_fm, mask_gm ] = meshgrid( mask_fs(1:end-1), mask_gs(1:end-1) );
    intensity2D = zeros(waferNf-1,waferNg-1,length(sourceData.Value));
    sourceX = sourceData.X;
    sourceY = sourceData.Y;
    sourceV = sourceData.Value;

    focus = SimulationRange(iFocus);
    dfmdg = (mask_fs(2)-mask_fs(1))*(mask_gs(2)-mask_gs(1));
    
    [ source_theta, source_rho] = cart2pol(sourceX,  sourceY);
    [ PolarizedX, PolarizedY ] = source.Calc_PolarizationMap( source_theta, source_rho );
    new_spectrum = spectrum(1:end-1,1:end-1);
    mask_fg2m = mask_fm.^2 + mask_gm.^2;
    sourceXY2 =  sourceX.^2+ sourceY.^2;
    
    for j = 1:length(sourceData.Value)  %%parfor
        obliqueRaysMatrix = 1;
        ExyzCalculateNumber_2D = 1;
        if strcmpi( numerics.ImageCalculationMode,'vector')
            ExyzCalculateNumber_2D = 3;
        elseif strcmpi( numerics.ImageCalculationMode,'scalar')
            ExyzCalculateNumber_2D = 1;
        end
        rho2 = mask_fg2m +2*(sourceX(j)*mask_fm + sourceY(j)*mask_gm) + sourceXY2(j);
        
        validPupil = find(rho2<=1);
        f_calc = mask_fm(validPupil) + sourceX(j);
        g_calc = mask_gm(validPupil) + sourceY(j);
        
        [theta_calc, rho_calc] = cart2pol(f_calc,g_calc);
        fgSquare = rho_calc.^2 ;
        
        obliquityFactor = sqrt(sqrt((1-(M^2*NA^2).*fgSquare)./ (1-((NA/indexImage)^2).*fgSquare)));
        aberration = projector.CalculateAberrationFast( rho_calc,theta_calc, Orientation);
        pupilFilter = projector.CalculatePupilFilter( rho_calc, theta_calc );%光瞳函数
        tempFocus = exp(-1j*2*pi/wavelength.*sqrt(indexK^2-NA*NA.*fgSquare).*focus);
        SpectrumCalc = new_spectrum(validPupil);
        
        TempHAber = SpectrumCalc.*obliquityFactor.*exp(1j*2*pi.*aberration).*pupilFilter.* tempFocus;%传递函数
        
        if strcmpi( numerics.ImageCalculationMode,'vector')
            obliqueRaysMatrix = zeros(length(fgSquare),ExyzCalculateNumber_2D);
            [M0xx,M0yx,M0xy,M0yy,M0xz,M0yz] = CalculateCharacteristicMatrix(f_calc,g_calc,fgSquare,NA,indexImage);
            obliqueRaysMatrix(:,1) = PolarizedX(j).*M0xx + PolarizedY(j).*M0yx;
            obliqueRaysMatrix(:,2) = PolarizedX(j).*M0xy + PolarizedY(j).*M0yy;
            obliqueRaysMatrix(:,3) = PolarizedX(j).*M0xz + PolarizedY(j).*M0yz;
        end
        
        rho2(:) = 0;
        intensityTemp = zeros(waferNf-1,waferNg-1);
        
        
        for iEM = 1:ExyzCalculateNumber_2D
            rho2(validPupil) = TempHAber.*obliqueRaysMatrix(:,iEM);
            
            if waferNf==maskNf && waferNg==maskNg
                ExyzFrequency = rho2;
            else
                ExyzFrequency = zeros(waferNf-1,waferNg-1);
                if waferNf > maskNf
                    rangeWaferNf = (waferNf-maskNf+2)/2 : (waferNf+maskNf-2)/2;
                    rangeMaskNf = 1:maskNf-1;
                else
                    rangeWaferNf =1:waferNf-1;
                    rangeMaskNf =  (maskNf-waferNf+2)/2 : (waferNf+maskNf-2)/2;
                end
                
                if waferNg > maskNg
                    rangeWaferNg =  (waferNg-maskNg+2)/2 : (waferNg+maskNg-2)/2;
                    rangeMaskNg = 1:maskNg-1;
                else
                    rangeWaferNg = 1:waferNg-1;
                    rangeMaskNg = (maskNg-waferNg+2)/2 : (waferNg+maskNg-2)/2;
                end
                ExyzFrequency(rangeWaferNf,rangeWaferNg) = rho2(rangeMaskNf,rangeMaskNg) ;
            end

            Exyz_Partial = fft2( ExyzFrequency );   %后续优化目标
            intensityTemp = intensityTemp + (real(Exyz_Partial).^2 + imag(Exyz_Partial).^2);

        end
        
        intensity2D(:,:,j) = intensityTemp;
    end
    intensity2D =  reshape(sourceV,1,1,[]).*intensity2D;
    intensity2D = dfmdg^2*fftshift(sum(intensity2D,3));
    
    
    intensity2D(:,waferNf) = intensity2D(:,1);
    intensity2D(waferNg,:) = intensity2D(1,:);
    intensity2D = real(rot90(intensity2D,2));
    Intensity(:,:,iFocus) = indexK*intensity2D/weight;
end
ImageX = linspace(-mask.Period_X/2,mask.Period_X/2,waferNf);
ImageY = linspace(-mask.Period_Y/2,mask.Period_Y/2,waferNg);
ImageZ = receipe.FocusRange;

lithoImage.ImageType = '2d';
lithoImage.SimulationType = 'aerial';
lithoImage.Intensity = Intensity;
lithoImage.ImageX =ImageX;
lithoImage.ImageY = ImageY;
lithoImage.ImageZ = ImageZ;

end