classdef FilmStack
    %% 光刻胶与基底
    % 支持光刻胶，顶层和底层光刻胶
    properties
        Layers
%         IndexResistReal = 1.7135;
%         IndexResistImag = 0.016894;
%         IndexSubstrateReal = 0.883143;
%         IndexSubstrateImag = 2.777792;
% 
%         ThinknessResist = 300;%Unit:nm

    end
    
    methods
        function fs =  FilmStack()
            l = Layer();
            l.Thickness = 400;
            l.Type = 'resist';
            l.Material = Resist.AR165J ;
            fs.Layers{1} = l;
            
            l.Thickness = inf;
            l.Type = 'substrate';
            l.Material = Material.Silicon;
            fs.Layers{2} = l;
            
        end
        
        function fs = AddLayer(fs,layer)%TODO:
            % 从上层往下层插入
            
        end
        
        function fs = RemoveLayer(fs,layerNo)%TODO: 待完成
            % 删除单层
            layerLength = length(fs.Layers);
            if layerNo >layerLength || layerNo < 1 
                error('Wrong layer Number !');
            end
            
            tempLayers = fs.Layers;
            fs.Layers = cell(1,layerLength-1);
            for iLayer = 1:length(tempLayers)
                if iLayer < layerNo
                    iLayer2 = iLayer;
                else
                    iLayer2 = iLayer + 1;
                end
                fs.Layers{iLayer} = tempLayers{iLayer2};
            end
        end
        
        function number = GetNumberOfLayers(fs)
            number = length(fs.Layers);
        end
        
        function layers = GetTARCLayers(fs)
            for i = 1:length(fs.Layers)
                if strcmpi(fs.Layers{i}.Type,'resist')
                    break;
                end
            end
            
            layers = cell(1,i-1);
            if i==1
                layers = [];
            elseif i<length(fs.Layers)
                for j = 1:i-1
                    layers{j} = fs.Layers{j};
                    layers{j}.IndexComplex = fs.Layers{j}.RefractiveIndex;
                end
            else
                error('No Resist!');
            end
        end
        
        function layer = GetResistLayer(fs)
            for i = 1:length(fs.Layers)
                if strcmpi(fs.Layers{i}.Type,'resist')
                    layer = fs.Layers{i};
                    return;
                end
            end
            error('No Resists Layer');
        end
        
        function layers = GetBARCLayers(fs)
            iResist = [];
            iSubstrate = [];
            for i = 1:length(fs.Layers)
                if strcmpi(fs.Layers{i}.Type,'resist')
                    iResist = i;
                elseif strcmpi(fs.Layers{i}.Type,'substrate')
                    iSubstrate = i;
                end
            end
            
            layers = cell(1,iSubstrate-iResist-1);
            if isempty(iResist)
                error('No Resist!');
            elseif isempty(iSubstrate)
                error('No Substrate!');
            elseif iSubstrate < iResist+1.5
               layers = [];
            elseif iResist < iSubstrate
                sulLayersList = iResist+1:iSubstrate-1;
                for i=1:length(layers)
                    layers{i} = fs.Layers{sulLayersList(i)};
                    layers{i}.IndexComplex = fs.Layers{sulLayersList(i)}.RefractiveIndex;
                end
            else
                error('Wrong filmstack setting');
            end
        end
        
        function layers = GetFullARCLayers(fs)
            layers = cell(1,length(fs.Layers)-1);
            for i = 1:length(fs.Layers)-1
                layers{i} = fs.Layers{i};
                 layers{i}.IndexComplex = fs.Layers{i}.RefractiveIndex;
            end
        end
        
        function index = GetResistIndex(fs)
            for i = 1:length(fs.Layers)
                if i == length(fs.Layers)
                    error('Wrong resist setting');
                elseif strcmpi(fs.Layers{i}.Type,'resist')
                    resistLayer = fs.Layers{i};
                    break;
                end
            end
            index = resistLayer.RefractiveIndex;
        end
        
        function index = GetSubstrateIndex(fs)
            substrateLayer = fs.Layers{end};
            if ~strcmpi( substrateLayer.Type,'substrate') 
                error('Wrong substrate setting');
            end
            index = substrateLayer.RefractiveIndex;
        end
        
        function thickness = GetResistThickness(fs)
            for i = 1:length(fs.Layers)
                if i == length(fs.Layers)
                    error('Wrong resist setting');
                elseif strcmpi(fs.Layers{i}.Type,'resist')
                    resistLayer = fs.Layers{i};
                    break;
                end
            end
            thickness = resistLayer.Thickness;
        end
    end
end