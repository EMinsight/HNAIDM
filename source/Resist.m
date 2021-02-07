classdef Resist < Material
   properties
       General_Parameter

       Unexposed_n
       Exposed_n
       Type = 'novolac' ; % 'novolac' / 'CAR'
       PEB_A = 0;     % unit:1/um  ��ѧ�Ŵ�A=0
       PEB_B;         % unit:1/um  JSR AR165J 1.1 ͨ������k����
       PEB_C = 0.05;  % unit:1/um cm2/mJ JSR AR165J 0.0434
       PEB_LnAr % ��ɢϵ�� nm2/s
       PEB_Ea  % ��� kcal/mol
       
       Development_Rmax;% �����Ӱ����
       Development_Rmin;% ��С��Ӱ����
       Development_Mth; % PACŨ����ֵ
       Development_Rresin;
       Development_n; % 
       Development_l; % inhibition reaction order
       Development_Rrs; %���淴Ӧ�������
       Development_Din; %��������ЧӦ���
       Tone = 'Positive' % ��̽����� 'Positive' / 'Negative'
   end
   
   methods 
        %���캯��
       function r = Resist()
           
       end
       
       function value = get.PEB_B(r)
          value = 4*pi*r.k/(r.WaveLength/1e3); %
       end
       
       function r = set.PEB_B(r, value)
           r.k = value * r.WaveLength/1e3 / (4*pi);
       end
       
   end
    
   methods (Static)
       function r = AR165J()
           r = Resist;
           r.Name = 'AR165J';
           r.n = 1.7135;
           r.k = 0.016894;
           r.WaveLength = 193;
           
           r.PEB_A = 0;
           r.PEB_C = 0.0434;
           
           r.PEB_LnAr = 27.12; % nm2/s
           r.PEB_Ea = 19.65; % ��� kcal/mol
           
           r.Development_Rmax = 100;
           r.Development_Rmin = 0.1;
           r.Development_Rresin = 3.5;
           r.Development_Mth = -100;
           r.Development_n = 3;
           r.Development_Rrs = 0.6;
           r.Development_Din = 120;
        end
   end
    
end