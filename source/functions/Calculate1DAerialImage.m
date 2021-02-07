function aerailLithoImage = Calculate1DAerialImage(sr, mk, po, aerailLithoImage, rp, numerics)
% �Ż��汾
maskNf = numerics.SampleNumber_Mask_X;  % 81
waferNf = numerics.SampleNumber_Wafer_X; % 81

sr.PntNum = numerics.SampleNumber_Source;
[ sourceData, ~ ] =  sr.Calc_SourceSimple();

weight = sum( sourceData.Value );
wavelength = sr.Wavelength; % 193
NA = po.NA;
M = po.Reduction;
if strcmpi(po.LensType,'Immersion')
    indexImage = po.Index_ImmersionLiquid; % Default = 1.44
elseif strcmpi(po.LensType,'Dry')
     indexImage = 1;  
     if NA>=1
         error('Wrong NA!');
     end
else
    error('Unsupported Lens Type')
end

% ��ģ��ת
if mk.Orientation == 0
    % �ޱ仯
elseif abs( mk.Orientation - pi/2)<eps
    tenpY = sourceData.Y;
    sourceData.Y = -sourceData.X;
    sourceData.X = tenpY;
else
    error('Not supported orientation angle')
end

% ��ģƵ��
mk.Nf = maskNf;
[ spectrum, mask_fs, ~ ] = mk.CalculateMaskSpectrum(po, sr);

%
SimulationRange = rp.FocusRange - rp.Focus ;

%% ��ʼ����
Intensity = zeros(length(SimulationRange),waferNf-1);%%%%%%%%

ExyzCalculateNumber = 1;

fm_s = repmat(sourceData.X,1,maskNf-1);
gm_s = repmat(sourceData.Y,1,maskNf-1);

fm_sm = fm_s + mask_fs(1,1:end-1); % NOTE:   %%��ҪMATLAB 2017���ϰ汾
gm_sm = gm_s;  % 

rho2 = fm_sm.^2 + gm_sm.^2;

validPupil = (rho2<=1);     % ͨ����һ���Ŀ׾����޶���Ч�׾�
f_calc = fm_sm(validPupil); 
g_calc = gm_sm(validPupil);

[theta_calc, rho_calc] = cart2pol(f_calc,g_calc);
fgSquare = rho_calc.^2 ;

% ��б����
obliquityFactor = sqrt(sqrt((1-(M^2*NA^2).*fgSquare)./(1-(NA/indexImage)^2.*fgSquare))); 
% ���㲨���
if mk.Orientation == 0
    % �ޱ仯
    aberration = po.CalculateAberrationFast(rho_calc,theta_calc, 0);
elseif abs( mk.Orientation - pi/2)<eps
    aberration = po.CalculateAberrationFast(rho_calc,theta_calc, pi/2);%�Ѿ��޸�
else
    error('Not supported orientation angle')
end

SpectrumCalc = ones(length(sourceData.Value),1)*spectrum(1:end-1); % TODO: ��Ҫ�Ż�
SpectrumCalc = SpectrumCalc(validPupil); % ��Ч�׾���Χ�ڵ���ģƵ��

TempH0Aber = SpectrumCalc.* obliquityFactor.*exp(1j*2*pi*aberration); 
% ��ģƵ�ס����������б���ӡ�����ȱ�ٹ�ͫ����������ģ����Դ����������ͶӰ�ﾵ��


obliqueRaysMatrix = ones(length(rho_calc),ExyzCalculateNumber); % ��б���߾���
if strcmp( numerics.ImageCalculationMode,'vector')
    ExyzCalculateNumber = 3;
    
    f_calc_s = fm_s(validPupil);
    g_calc_s = gm_s(validPupil);
    [ nts, nrs ] = cart2pol(f_calc_s,g_calc_s);
    
    %ƫ��ֽ�
    [ PolarizedX, PolarizedY ] = sr.Calc_PolarizationMap( nts, nrs );
    
    % ����ƫ�������M
    [M0xx,M0yx,M0xy,M0yy,M0xz,M0yz] = CalculateCharacteristicMatrix(f_calc,g_calc,fgSquare,NA,indexImage);
        
    obliqueRaysMatrix = ones(length(rho_calc),ExyzCalculateNumber);  % ����ʸ��ģ������������б���߾���
    obliqueRaysMatrix(:,1) = PolarizedX.*M0xx + PolarizedY.*M0yx;
    obliqueRaysMatrix(:,2) = PolarizedX.*M0xy + PolarizedY.*M0yy;
    obliqueRaysMatrix(:,3) = PolarizedX.*M0xz + PolarizedY.*M0yz;
    
end


if ~strcmp(po.PupilFilter.Type,'none') % ��ʱ��ͫ�������ڵ�ͨ�˲����ԣ���������pupilFilterGauss
    filter = po.PupilFilter.Type;
    parameter = po.PupilFilter.Parameter;
    f_pupil = fm_sm(validPupil);
    g_pupil = gm_sm(validPupil);
    pupilFilterData = feval(filter,parameter,f_pupil,g_pupil); % ��ͫ���������������������Ϊ�������������Ƶ�������ṩ����
    TempH0Aber = pupilFilterData.*TempH0Aber; % ������Ƶ�׼���
    
end

TempFocus = -1j*2*pi/wavelength.*sqrt(indexImage^2 - NA*NA.*fgSquare); % ��������Ǳ�ʾ�뽹����

for iFocus = 1:length(SimulationRange)
    intensity = zeros(1,waferNf-1);

    
    tempF = exp(TempFocus.*SimulationRange(iFocus));
    for iEM = 1:ExyzCalculateNumber
        rho2(:) = 0;
        rho2(validPupil) = obliqueRaysMatrix(:,iEM).*TempH0Aber .*tempF;
        if waferNf == maskNf
            ExyzFrequency = rho2;
        elseif waferNf >maskNf
            ExyzFrequency = [zeros( length(sourceData.Value), (waferNf-maskNf)/2),rho2,zeros( length(sourceData.Value), (waferNf-maskNf)/2)];
        else
            ExyzFrequency = rho2;
            ExyzFrequency(:, 1:(maskNf-waferNf)/2) = [ ];
            ExyzFrequency(:, (end-(maskNf-waferNf)/2+1):end) = [ ];
        end

        Exyz = fft(ExyzFrequency,[],2); 
        IntensityCon = real(Exyz).^2 + imag(Exyz).^2;
        
        IntensityTemp = indexImage * (mask_fs(2) - mask_fs(1)).^2 * fftshift(sourceData.Value' * IntensityCon);
        intensity = intensity + IntensityTemp;

    end

    Intensity(iFocus,:) = intensity./weight;

end
ImageX = linspace(-mk.Period_X/2,mk.Period_X/2,waferNf);
ImageY = 0;
ImageZ = rp.FocusRange;

aerailLithoImage.ImageType = '1d';
aerailLithoImage.SimulationType = 'aerial';
aerailLithoImage.Intensity = fliplr([Intensity, Intensity(:,1)]);
aerailLithoImage.ImageX = ImageX;
aerailLithoImage.ImageY = ImageY;
aerailLithoImage.ImageZ = ImageZ;

end