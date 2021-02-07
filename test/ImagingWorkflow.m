%% Abbe imaging test
%% 1DMask
clear
close all;

scanner = ImagingModel;
scanner.Source.Shape = 'dipolecirc';
scanner.Source.SigmaOut = 0.9;
scanner.Source.SigmaIn = 0.7;
scanner.Source.OpenAngle = pi/6;
scanner.Source.AdvancedParametersEnable = 1;
% lm.Projector.Aberration_Zernike = rand(1,37)*0.02;
% scanner.Mask = Mask.CreateMask('line_space');

scanner.Numerics.ImageCalculationMode = 'vector';
scanner.Source.PolarizationType = 'y_pol';

scanner.Mask = Mask.CreateMask('space',200,600);


% sr = scanner.Source;
% po = scanner.Projector;
% mk = scanner.Mask;
% rp = scanner.Receipe;
% nm = scanner.Numerics;
% aliN = CalculateNormalImage(sr,mk,po,ali,rp,nm);

figure;
ImagingView.PlotSource(scanner.Source);
figure;
ImagingView.PlotMask(scanner.Mask);
figure;
ImagingView.PlotMaskSpectrum(scanner.Mask,scanner.Source,scanner.Projector);

%% Abbe Imaging
% 1D mask
scanner.Numerics.Normalization_Intensity = 0;
ali = scanner.CalculateAerialImage();
figure;
ImagingView.PlotAerialImage(ali);

scanner.Numerics.Normalization_Intensity = 1;
ali = scanner.CalculateAerialImage();
figure;
ImagingView.PlotAerialImage(ali);

% 2D Mask
scanner.Mask = Mask.CreateMask('line_space');
aliA = scanner.CalculateAerialImage();
figure;
ImagingView.PlotAerialImage(aliA);

%% Hopkins Imaging
scanner.Numerics.ImageCalculationMethod = 'Hopkins';
% scanner.Numerics.Hopkins_SettingType = 'order';
% scanner.Numerics.Hopkins_Order = 30;
aliH = scanner.CalculateAerialImage();
figure;
ImagingView.PlotAerialImage(aliH);

figure;
pcolor(aliA.Intensity-aliH.Intensity);
shading flat;
colorbar;

%% ½ºÄÚÏñ
scanner.Numerics.ImageCalculationMethod = 'Abbe';
%1d
scanner.Mask = Mask.CreateMask('space',200,600);
scanner.Receipe.FocusRange = [];
rli = scanner.CalculateResistImage();
figure;
ImagingView.PlotResistImage(rli);

%2d
scanner.Mask = Mask.CreateMask('line_space');
scanner.Receipe.FocusRange = 0;
rli = scanner.CalculateResistImage();
figure;
ImagingView.PlotResistImage(rli);


%% metrology

