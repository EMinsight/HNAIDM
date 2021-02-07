classdef Layer
    properties
        Type   % 'arc' or 'resist' or 'substrate' 无大小写要求
        Thickness % unit:nm
        Material  %Type: Material or Resist
        RefractiveIndex % Complex index
    end
    
    properties(Hidden = true)
        IndexComplex;
    end
    
    methods
        function l = Layer(type,thickness,material)
            switch nargin
                case 3
                    
                    l.Type = type;
                    l.Thickness = thickness;
                    l.Material = material;
                otherwise
            end
        end
        
        function index = get.RefractiveIndex(l)
            index = l.Material.n + 1j*l.Material.k;
% %             index = l.Material.n;
        end
        
        %         function index = RefractiveIndex(l)
        %             index = l.Material.n + 1j*l.Material.k;
        %         end
        
    end
end