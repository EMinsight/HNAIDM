classdef ProjectionObjective
    %% properities
    properties
        Aberration_Zernike;
        Aberration_Polarization;
        PupilFilter;%'none', 'gaussian'
        Reduction = 0.25;
        NA = 1.35;% 0.65-0.93 for dry mode, 1.0-1.35 for immersion mode
        LensType = 'Immersion'% Immersion / Dry 
        Index_ImmersionLiquid = 1.44; % Air => 1.0, Water => 1.44 
    end
    
    %% functions
    methods
        function p = ProjectionObjective()
            p.PupilFilter.Type = 'none';
            p.PupilFilter.Parameters = 0;
            p.Aberration_Zernike = zeros(1,37);
        end
        
        function aberration = CalculateAberration(p, rho, theta,Orientation)
            
            %Rotate projection objective            
            if abs(Orientation)>eps
                Coefficients = RoateAngle(p.Aberration_Zernike,Orientation);
            else
                Coefficients = p.Aberration_Zernike;
            end
            
            % Calculate aberrations
            aberration= Coefficients(1) .* ones(size(theta))...
                +Coefficients(2) .* rho.*cos(theta) ...
                +Coefficients(3) .* rho.*sin(theta)...
                +Coefficients(4) .* ( 2*rho.*rho-1)...
                +Coefficients(5) .* rho.^2.*cos(2*theta)...
                +Coefficients(6) .* rho.^2.*sin(2*theta)...
                +Coefficients(7) .* (3*rho.^3-2*rho).*cos(theta)...
                +Coefficients(8) .* (3*rho.^3-2*rho).*sin(theta)...
                +Coefficients(9) .* ( 6*rho.^4-6*rho.^2+1)...
                +Coefficients(10).* rho.^3.*cos(3*theta)...
                +Coefficients(11).* rho.^3.*sin(3*theta)...
                +Coefficients(12).* (4*rho.^4-3*rho.^2).*cos(2*theta)...
                +Coefficients(13).* (4*rho.^4-3*rho.^2).*sin(2*theta)...
                +Coefficients(14).* (10*rho.^5-12*rho.^3+3*rho).*cos(theta)...
                +Coefficients(15).* (10*rho.^5-12*rho.^3+3*rho).*sin(theta)...
                +Coefficients(16).* (20*rho.^6-30*rho.^4+12*rho.^2-1)...
                +Coefficients(17).* (rho.^4).*cos(4*theta)...
                +Coefficients(18).* (rho.^4).*sin(4*theta)...
                +Coefficients(19).* (5*rho.^5-4*rho.^3).*cos(3*theta)...
                +Coefficients(20).* (5*rho.^5-4*rho.^3).*sin(3*theta)...
                +Coefficients(21).* (15*rho.^6-20*rho.^4+6*rho.^2).*cos(2*theta)...
                +Coefficients(22).* (15*rho.^6-20*rho.^4+6*rho.^2).*sin(2*theta)...
                +Coefficients(23).* (35*rho.^7-60*rho.^5+30*rho.^3-4*rho).*cos(theta)...
                +Coefficients(24).* (35*rho.^7-60*rho.^5+30*rho.^3-4*rho).*sin(theta)...
                +Coefficients(25).* (70*rho.^8-140*rho.^6+90*rho.^4-20*rho.^2+1)...
                +Coefficients(26).* (rho.^5).*cos(5*theta)...
                +Coefficients(27).* (rho.^5).*sin(5*theta)...
                +Coefficients(28).* (6*rho.^6-5*rho.^4).*cos(4*theta)...
                +Coefficients(29).* (6*rho.^6-5*rho.^4).*sin(4*theta)...
                +Coefficients(30).* (21*rho.^7-30*rho.^5+10*rho.^3).*cos(3*theta)...
                +Coefficients(31).* (21*rho.^7-30*rho.^5+10*rho.^3).*sin(3*theta)...
                +Coefficients(32).* (56*rho.^8-105*rho.^6+60*rho.^4-10*rho.^2).*cos(2*theta)...
                +Coefficients(33).* (56*rho.^8-105*rho.^6+60*rho.^4-10*rho.^2).*sin(2*theta)...
                +Coefficients(34).* (126*rho.^9-280*rho.^7+210*rho.^5-60*rho.^3+5*rho).*cos(theta)...
                +Coefficients(35).* (126*rho.^9-280*rho.^7+210*rho.^5-60*rho.^3+5*rho).*sin(theta)...
                +Coefficients(36).* ( 252*rho.^10-630*rho.^8+560*rho.^6-210*rho.^4+30*rho.^2-1 ) ...
                +Coefficients(37).* (924*rho.^12-2772*rho.^10+3150*rho.^8-1680*rho.^6+420*rho.^4-42*rho.^2+1);
            
        end
        
        function aberration = CalculateAberrationFast(p, rho, theta,Orientation)
            % Perforamance optimized version
            % radial value
            r2 = rho.^2;
            r3 = rho.* r2;
            r4 = r2.^2;
            r5 = rho.* r4;
            r6 = r3.^2;
            r7 = rho.* r6;
            r8 = r4.^2;
            r9 = rho.* r8;
            r10 = r5.^2;
            r12 = r6.^2;
            
            % azimuthal value
            ct = cos(theta);
            st = sin(theta);
            
            c2t = ct.*ct-st.*st;
            s2t = 2*ct.*st;
            
            c3t = c2t.*ct-s2t.*st;
            s3t = s2t.*ct+c2t.*st;
            
            c4t = c3t.*ct-s3t.*st;
            s4t = s3t.*ct+c3t.*st;
            
            c5t = c4t.*ct-s4t.*st;
            s5t = s4t.*ct+c4t.*st;
            
            if abs(Orientation)>eps
                Coefficients = ProjectionObjective.RoateAngle(p.Aberration_Zernike,Orientation);
            else
                Coefficients = p.Aberration_Zernike;
            end
            
            %calcualte aberration distribution
            aberration = Coefficients(1);
            aberration = aberration + Coefficients(2) .* rho.*ct;
            aberration = aberration + Coefficients(3) .* rho.*st;
            aberration = aberration + Coefficients(4) .* ( 2*r2-1);
            aberration = aberration + Coefficients(5) .* r2.*c2t;
            aberration = aberration + Coefficients(6) .* r2.*s2t;
            aberration = aberration + Coefficients(7) .* (3*r3-2*rho).*ct;
            aberration = aberration + Coefficients(8) .* (3*r3-2*rho).*st;
            aberration = aberration + Coefficients(9) .* ( 6*r4-6*r2+1);
            aberration = aberration + Coefficients(10).* r3.*c3t;
            aberration = aberration + Coefficients(11).* r3.*s3t;
            aberration = aberration + Coefficients(12).* (4*r4-3*r2).*c2t;
            aberration = aberration + Coefficients(13).* (4*r4-3*r2).*s2t;
            aberration = aberration + Coefficients(14).* (10*r5-12*r3+3*rho).*ct;
            aberration = aberration + Coefficients(15).* (10*r5-12*r3+3*rho).*st;
            aberration = aberration + Coefficients(16).* (20*r6-30*r4+12*r2-1);
            aberration = aberration + Coefficients(17).* r4.*c4t;
            aberration = aberration + Coefficients(18).* r4.*s4t;
            aberration = aberration + Coefficients(19).* (5*r5-4*r3).*c3t;
            aberration = aberration + Coefficients(20).* (5*r5-4*r3).*s3t;
            aberration = aberration + Coefficients(21).* (15*r6-20*r4+6*r2).*c2t;
            aberration = aberration + Coefficients(22).* (15*r6-20*r4+6*r2).*s2t;
            aberration = aberration + Coefficients(23).* (35*r7-60*r5+30*r3-4*rho).*ct;
            aberration = aberration + Coefficients(24).* (35*r7-60*r5+30*r3-4*rho).*st;
            aberration = aberration + Coefficients(25).* (70*r8-140*r6+90*r4-20*r2+1);
            aberration = aberration + Coefficients(26).* (r5).*c5t;
            aberration = aberration + Coefficients(27).* (r5).*s5t;
            aberration = aberration + Coefficients(28).* (6*r6-5*r4).*c4t;
            aberration = aberration + Coefficients(29).* (6*r6-5*r4).*s4t;
            aberration = aberration + Coefficients(30).* (21*r7-30*r5+10*r3).*c3t;
            aberration = aberration + Coefficients(31).* (21*r7-30*r5+10*r3).*s3t;
            aberration = aberration + Coefficients(32).* (56*r8-105*r6+60*r4-10*r2).*c2t;
            aberration = aberration + Coefficients(33).* (56*r8-105*r6+60*r4-10*r2).*s2t;
            aberration = aberration + Coefficients(34).* (126*r9-280*r7+210*r5-60*r3+5*rho).*ct;
            aberration = aberration + Coefficients(35).* (126*r9-280*r7+210*r5-60*r3+5*rho).*st;
            aberration = aberration + Coefficients(36).* (252*r10-630*r8+560*r6-210*r4+30*r2-1 );
            aberration = aberration + Coefficients(37).* (924*r12-2772*r10+3150*r8-1680*r6+420*r4-42*r2+1);
        end
            
        function pupilFilter = CalculatePupilFilter(p,rho,theta)
            parameters = p.PupilFilter;
            if strcmpi(parameters.Type,'gaussian')
                filter = p.PupilFilter.FilterType;
                para = parameters.FilterPara;
                pupilFilter = feval(filter,para,rho,theta);   
            elseif strcmpi(p.PupilFilter.Type,'none')
                pupilFilter = 1;
            else
                pupilFilter = 1;
            end
        end
    end
    
    
    methods (Static)
        function c1 = RoateAngle(c0,theta) %Rotate zernike aberration
            %对于不在坐标轴的一维掩模，其像差与在坐标轴上的一维掩模的像差关系
            %不含cos、sin项的无变化
            %含有cos项的等于本身乘以cos加上对应的含sin项的乘以sin
            %含有sin项的等于本身乘以sin加上对应的含cos项的乘以cos
            mm=[0 1 1 0 2 2 1 1 0 3 3 2 2 1 1 0 4 4 3 3 2 2 1 1 0 5 5 4 4 3 3 2 2 1 1 0 0];%m
            tt=[0 1 -1 0 1 -1 1 -1 0 1 -1 1 -1 1 -1 0 1 -1 1 -1 1 -1 1 -1 0 1 -1 1 -1 1 -1 1 -1 1 -1 0 0];
            pp=[1 3 2 4 6 5 8 7 9 11 10 13 12 15 14 16 18 17 20 19 22 21 24 23 25 27 26 29 28 31 30 33 32 35 34 36 37];
            
            % theta=theta*pi/180; %角度转为弧度 -弧度制输入
            c1=zeros(1,37);
            
            for ii=1:37
                c1(ii) = c1(ii) + c0(ii)*cos( mm(ii)* theta );
                c1( pp(ii) ) = c1( pp(ii) ) - tt(ii)* c0(ii)*sin( mm(ii)* theta );
            end
        end
    end
end
