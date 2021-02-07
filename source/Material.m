classdef Material
    properties
        Name
        n
        k
        WaveLength
    end
    methods
        function m = Material(name,indexn,indexlk,wavelength)
            switch nargin
                case 4
                    m.Name = name;
                    m.n = indexn;
                    m.k = indexlk;
                    m.WaveLength = wavelength;
                otherwise
            end
        end
    end
    
    methods (Static)
        function m = Silicon()
            m = Material('Si',0.883143,2.777792,193);
        end
        
        function m = SiliconDioxide()
            m = Material;
            m.Name = 'SiO2';
            m.n = 1.563117;
            m.k = 0;
            m.WaveLength = 193;
        end
        
        
        function m = ARCTop()
            m = Material;
            m.Name = 'AZ Aquatar ArF Top ARC';
            m.n = 1.513282;
            m.k = 0.004249;
            m.WaveLength = 193;
        end
        
        function m = ARCBottom()
            m = Material;
            m.Name = 'Brewer DUV 42C';
            m.n = 1.48;
            m.k = 0.41;
            m.WaveLength = 193;
        end
    end
end