classdef Numerics
    %%%%%%%%%%%%%%%%%
    % Universal parameter for model
    %
    %%%%%%%%%%%%%%%%%
    properties
        SampleNumber_Source = 101; 
        
        SampleNumber_Mask_X = 81;
        SampleNumber_Mask_Y = 81;
        SampleNumber_Mask_RCWA_X = 11; % Sample number for RCWA Mask
        SampleNumber_Mask_RCWA_Y = 11; %
        
        Sample_Calc_Mode = 'mask'; % sampling number for mask calculatuon => mask / auto
        
        SampleNumber_Wafer_X = 81;
        SampleNumber_Wafer_Y = 81;
        SampleNumber_Wafer_Z = 41;
        
        % simulation settings
        SimulationRange_Aerial = 0;
        SimulationRange_Resist = [];
        % Normailization
        Normalization_Intensity = 1;
        
        ImageCalculationMode = 'vector'; % Imaging model => 'vector'  'sclar' 
        ImageCalculationMethod = 'abbe'; % Calculation model =>¡®abbe¡¯ ¡®hopkins¡¯
        Hopkins_SettingType = 'order'; % Truncation method for hopkins model => 'order' 'Threshold'
        Hopkins_Order = 50;
        Hopkins_Threshold = 0.95 % (0,1)
        
        ResisitModel = 'physical'; % Resist model => physical / lumped / measured        
        DevelopmentModel = 'Mack'% Development model => Mack / Enhanced Mack / Notch        
        
    end
    
    properties (Constant) 
        % Physical constant
        Mu0 = 4*pi*10-7; %Permeability => T¡¤m/A
        Epsilon0 = 8.854187817*10-12; %Permittivity => F/m
        R = 0.0019876 % Gas constant => kcal/(K*mol)
        K0 = 273.15 %Absolute zero => K
    end
    
    methods
        function n =  Numerics()
            
        end
    end
    
end
