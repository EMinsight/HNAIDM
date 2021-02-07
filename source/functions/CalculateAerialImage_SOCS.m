function intensity = CalculateAerialImage_SOCS(mask,TCCMatrix_SOCS,source,projector,numerics)

%% Get Image
maskNf = numerics.SampleNumber_Mask_X;  % 默认81
maskNg = numerics.SampleNumber_Mask_Y;  % 默认81

waferNf = numerics.SampleNumber_Wafer_X;  % 默认81
waferNg = numerics.SampleNumber_Wafer_Y;  % 默认81

[ spectrum, f, g ] = mask.CalculateMaskSpectrum(projector ,source); 


if waferNf==maskNf && waferNg==maskNg
     spectrumEx = spectrum(1:end-1,1:end-1);
else
    
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
    
    spectrumEx = zeros(waferNf-1,waferNg-1);
    spectrumEx(rangeWaferNf,rangeWaferNg) = spectrum(rangeMaskNf,rangeMaskNg) ;
    
end

temp = TCCMatrix_SOCS.*fftshift(spectrumEx);
Etemp = (f(2)-f(1))*(g(2)-g(1))*fft2(temp);
% Esquare = real(Etemp).^2 + imag(Etemp).^2;
Esquare = abs(Etemp).^2;

intensity =  sum(Esquare,3);
intensity = fftshift(gather(intensity));
intensity(:,waferNf) = intensity(:,1);
intensity(waferNg,:) = intensity(1,:);

intensity = rot90(intensity,2);

end