classdef Receipe
    properties
        % Exposure settings
        Dose = 30; % Unit: mJ/cm2
        
        % focus setting
        % use the top surface of FilmStack as reference
        % negative -> wafer 
        % positive -> projection objective
        Focus = 0; % Unit: nm
        FocusRange = 0;
        FocusReference = 0; % Reference offset
        
        % PEB²ÎÊý
        PEB_Temperature = 120; % centigrade ¡æ
        PEB_Time = 90 % PEB time -> second
        
        % Developing parameters

    end
    
end