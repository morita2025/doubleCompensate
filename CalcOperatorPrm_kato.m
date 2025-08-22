classdef CalcOperatorPrm_kato < handle
   properties (Constant ,Access = private)


        %d(1)-%d(9)
        d = [0.12; 0.07; 0.03; 0.03; 0.02; 0.01; 0.01; 0.02;];
        dx13 = CalcOperatorPrm_kato.d(3,1) + CalcOperatorPrm_kato.d(6,1) + CalcOperatorPrm_kato.d(7,1);
        dx2  = CalcOperatorPrm_kato.d(8,1);
        dxThermalCond = (0.5*CalcOperatorPrm_kato.d(3,1) + CalcOperatorPrm_kato.d(7,1) + 0.5*CalcOperatorPrm_kato.d(8,1) )*ones(3,1);


        %S(1)-S(9)
        S = [2.6e-3; 7e-4; 9.8e-3; 9e-4; pi*9e-6; pi*3e-4; 1.4e-3; 2.8e-4; pi*1.2e-4;];

        
        %aliminium 密度  %熱伝導率 %比熱 %放射率
        aluminum = struct("density",2500,"thermalCond",150,"specificHeat",250,"emissivity",0.2);

        
        %water密度  %熱伝導率 λ   %比熱  %放射率 ε
        water = struct("density",998.2,"thermalCond",0.602,"specificHeat",4182*0.52,"emissivity",0.93);
        
        
        %quantity
        quantity = struct("boltzmannConst",5.67*(10^(-8)),"absoluteTemperature",273.15);


        %peltier
        peltier =  struct("seebeck",0.2,"thermalConductance",0.63,...
                "resistance",2,"absoluteTemperature",273.15);
        

        end
        



   properties (SetAccess = immutable) 
        kp;
        ki;
        A_a;
        A_w;


        G_a;
        G_w;
        G_aw;
        ma_ca13; 
        ma_ca2;
        mw_cw13;
        mw_cw2;

        T_omega;
        T_0;

        interferrenceConstPrm;
        MOperatorConstPrm;
        peltierPrm;
        lowPassFilterTimePrm;
        tubeInterferrenceConstPrm;





   end

   properties (Access = public)
       settings;
   end






    methods

        function obj = CalcOperatorPrm_kato(options)
        
            arguments
                options.isRugekuttaMethodUse = false;
                options.isinterferenceRejectionUse = true;
                options.isFaltTolerantControlUse = false;
                options.isInterferrence = true;
                options.isD1Compensate =true;
                options.isD2Compensate=false;
                options.max_time = 0;


                %options.simulationTime {mustBeNumeric} = 600;
                options.outsideTemperature {mustBeNumeric} = 22;  %T_0
                options.heatTransferCoef_air  {mustBeNumeric} = 150; %α
                options.heatTransferCoef_water {mustBeNumeric} = 440;%α_ω 
   

                options.kp {mustBeNumeric} = [0; 0];
                options.ki {mustBeNumeric} = [0; 0];
                options.k {mustBeNumeric} = [1; 1];
                options.p {mustBeNumeric} = 1;
                options.p2 {mustBeNumeric} = 1;
                options.tau {mustBeNumeric} = 5;

                %current
                options.i_max = 1;
                options.i_min = -0.5;





            end
            obj.settings = options;
            obj.settings.initialWaterTemperature  = obj.settings.outsideTemperature;
            
           


            %パラメータ計算
            %A_nm--------------------------------------------------------------------------------
            A_a = zeros(3,4);
            m_a13 = obj.dx13 *( obj.S(3,1) - obj.S(5,1)) * obj.aluminum.density;
            m_a2 =  obj.dx2 *( obj.S(3,1) - obj.S(5,1)) * obj.aluminum.density;
            
            
            %Aa_11
            b1 = 2*obj.peltier.thermalConductance;
            b2 = obj.settings.heatTransferCoef_air* (2*obj.S(1,1) + 2*obj.S(2,1) + obj.S(3,1) -obj.S(5,1));
            b3 = obj.settings.heatTransferCoef_water*obj.S(6,1);
            b4 = obj.aluminum.thermalCond*obj.S(3,1)/obj.dxThermalCond(1);
            b5 = 4*( (obj.settings.outsideTemperature+obj.quantity.absoluteTemperature)^3 )*obj.aluminum.emissivity * obj.quantity.boltzmannConst*(2*obj.S(1,1) + 2*obj.S(2,1) + obj.S(3,1) -obj.S(5,1));
            c1 = m_a13 * obj.aluminum.specificHeat; 
            
            A_a(1,1) =(b1 + b2 + b3 + b4 + b5) / c1;
            
            %Aa_12
            b1 = 6*obj.aluminum.emissivity * obj.quantity.boltzmannConst * (( obj.settings.outsideTemperature+obj.quantity.absoluteTemperature)^2 );
            b2 = obj.S(1,1) + obj.S(2,1) + obj.S(3,1) - 2*obj.S(4,1) -obj.S(5,1);
            c1 = m_a13 * obj.aluminum.specificHeat; 
            
            A_a(1,2) =(b1 * b2 ) / c1;
            
            %Aa_13
            b1 = 4*obj.aluminum.emissivity * obj.quantity.boltzmannConst * ((obj.settings.outsideTemperature+obj.quantity.absoluteTemperature)^(1));
            b2 = obj.S(1,1) + obj.S(2,1) + obj.S(3,1) - 2*obj.S(4,1) -obj.S(5,1);
            c1 = m_a13 * obj.aluminum.specificHeat; 
            
            A_a(1,3) =(b1 * b2 ) / c1;
            
            %Aa_14
            b1 = obj.aluminum.emissivity * obj.quantity.boltzmannConst ;
            b2 = obj.S(1,1) + obj.S(2,1) + obj.S(3,1) - 2*obj.S(4,1) -obj.S(5,1);
            c1 = m_a13 * obj.aluminum.specificHeat; 
            
            A_a(1,4) =(b1 * b2 ) / c1;
            
            %Aa_21
            b1 = obj.settings.heatTransferCoef_air + 4*obj.aluminum.emissivity * obj.quantity.boltzmannConst * ((obj.settings.outsideTemperature+obj.quantity.absoluteTemperature)^(3));
            b2 = 2*obj.S(7,1) + 2*obj.S(8,1);
            b3 = obj.settings.heatTransferCoef_water * obj.S(9,1);
            b4 = (2 * obj.aluminum.thermalCond * obj.S(3,1)) / obj.dxThermalCond(2) ;
            c1 = m_a2 * obj.aluminum.specificHeat; 
            
            A_a(2,1) = (b1 * b2 + b3 + b4) / c1;
            
            %Aa_22
            b1 = 6*obj.aluminum.emissivity * obj.quantity.boltzmannConst * ((obj.settings.outsideTemperature+obj.quantity.absoluteTemperature)^(2));
            b2 = 2*obj.S(7,1) + 2*obj.S(8,1);
            c1 = m_a2 * obj.aluminum.specificHeat; 
            
            A_a(2,2) = (b1 * b2) / c1;
            
            
            %Aa_23
            b1 = 4*obj.aluminum.emissivity * obj.quantity.boltzmannConst * ((obj.settings.outsideTemperature+obj.quantity.absoluteTemperature)^(1));
            b2 = 2*obj.S(7,1) + 2*obj.S(8,1);
            c1 = m_a2 * obj.aluminum.specificHeat; 
            
            A_a(2,3) = (b1 * b2) / c1;
            
            %A_24
            b1 = obj.aluminum.emissivity * obj.quantity.boltzmannConst;
            b2 = 2*obj.S(7,1) + 2*obj.S(8,1);
            c1 = m_a2 * obj.aluminum.specificHeat; 
            
            A_a(2,4) = (b1 * b2) / c1;
            
            %Aa_31
            b1 = 2*obj.peltier.thermalConductance;
            b2 = obj.settings.heatTransferCoef_air * (2*obj.S(1,1) + 2*obj.S(2,1) + obj.S(3,1) -obj.S(5,1));
            b3 = obj.settings.heatTransferCoef_water*obj.S(6,1);
            b4 = obj.aluminum.thermalCond*obj.S(3,1)/obj.dxThermalCond(3);
            b5 = 4*((obj.settings.outsideTemperature+obj.quantity.absoluteTemperature)^(3))*obj.aluminum.emissivity * obj.quantity.boltzmannConst*(2*obj.S(1,1) + 2*obj.S(2,1) + obj.S(3,1) -obj.S(5,1));
            c1 = m_a13 * obj.aluminum.specificHeat; 
            
            A_a(3,1) =(b1 + b2 + b3 + b4 + b5) / c1;
            
            %Aa_32
            b1 = 6*obj.aluminum.emissivity * obj.quantity.boltzmannConst * ((obj.settings.outsideTemperature+obj.quantity.absoluteTemperature)^(2));
            b2 = obj.S(1,1) + obj.S(2,1) + obj.S(3,1) - 2*obj.S(4,1) -obj.S(5,1);
            c1 = m_a13 * obj.aluminum.specificHeat; 
            
            A_a(3,2) =(b1 * b2 ) / c1;
            
            %A_33
            b1 = 4*obj.aluminum.emissivity * obj.quantity.boltzmannConst * ((obj.settings.outsideTemperature+obj.quantity.absoluteTemperature)^(1));
            b2 = obj.S(1,1) + obj.S(2,1) + obj.S(3,1) - 2*obj.S(4,1) -obj.S(5,1);
            c1 = m_a13 * obj.aluminum.specificHeat; 
            
            A_a(3,3) =(b1 * b2 ) / c1;
            
            %Aa_34
            b1 = obj.aluminum.emissivity * obj.quantity.boltzmannConst ;
            b2 = obj.S(1,1) + obj.S(2,1) + obj.S(3,1) - 2*obj.S(4,1) -obj.S(5,1);
            c1 = m_a13 * obj.aluminum.specificHeat; 
            
            A_a(3,4) =(b1 * b2 ) / c1;
            

            %代入
            obj.A_a = A_a;


            %A_w_mn-----------------------------------------------------------------------------------
            A_w = zeros(3,1);
            m_omega13  = obj.dx13 * obj.S(5,1) * obj.water.density;
            m_omega2  = obj.dx2 * obj.S(5,1) * obj.water.density;
            
            
            %Aw_1
            b1 = (2*obj.water.thermalCond * obj.S(5,1)) / obj.dx13;
            b2 = obj.settings.heatTransferCoef_water * obj.S(6,1);
            c1 = m_omega13 * obj.water.specificHeat;
            
            A_w(1,1) = (b1 + b2) / c1;
            
            %Aw_2
            b1 = (obj.water.thermalCond * obj.S(5,1)) / obj.dx2;
            b2 = obj.settings.heatTransferCoef_water * obj.S(9,1);
            c1 = m_omega2 * obj.water.specificHeat;
            
            A_w(2,1) = (b1 + b2) / c1;
            
            %Aw_3
            b1 = (2*obj.water.thermalCond * obj.S(5,1)) / obj.dx13;
            b2 = obj.settings.heatTransferCoef_water * obj.S(6,1);
            c1 = m_omega13 * obj.water.specificHeat;
            
            A_w(3,1) = (b1 + b2) / c1;
            
            
            %代入
            obj.A_w = A_w;


            %G-----------------------------------------------------------------
            % G_a = zeros(3,1);
            % G_a(1,1) = obj.aluminum.thermalCond * obj.S(3,1) / obj.dx2;
            % G_a(2,1) = obj.aluminum.thermalCond * obj.S(3,1) / obj.dx13;%kaeteru itijiteki
            % G_a(3,1) = obj.aluminum.thermalCond * obj.S(3,1) / obj.dx2;
            % 
            % G_w = zeros(3,1);
            % G_w(1,1) = obj.water.thermalCond * obj.S(5,1) / obj.dx2;
            % G_w(2,1) = 2 * obj.water.thermalCond * obj.S(5,1) / obj.dx13;
            % G_w(3,1) = obj.water.thermalCond * obj.S(5,1) / obj.dx2;
            % 
            % G_aw = zeros(3,1);
            % G_aw(1,1) = obj.settings.heatTransferCoef_water * obj.S(6,1);
            % G_aw(2,1) = obj.settings.heatTransferCoef_water * obj.S(9,1);
            % G_aw(3,1) = obj.settings.heatTransferCoef_water * obj.S(6,1);
             
            
            %unko
            % obj.G_a = G_a;
            % obj.G_w = G_w;
            % obj.G_aw = G_aw;
            obj.ma_ca13 = m_a13 * obj.aluminum.specificHeat;
            obj.ma_ca2 = m_a2 * obj.aluminum.specificHeat;
            obj.mw_cw13 = m_omega13 * obj.water.specificHeat;
            obj.mw_cw2 = m_omega2* obj.water.specificHeat;    
            obj.T_omega = obj.settings.initialWaterTemperature;
            obj.T_0 = obj.settings.outsideTemperature +obj.quantity.absoluteTemperature;
            



            obj.kp = [obj.settings.kp(1,1); 0; obj.settings.kp(2,1)];
            obj.ki = [obj.settings.ki(1,1); 0; obj.settings.ki(2,1)];


            obj.interferrenceConstPrm = struct("y_a",obj.aluminum.thermalCond*obj.S(3,1) ./ obj.dxThermalCond,...
                                               "y_aw",obj.settings.heatTransferCoef_water.*[obj.S(6,1);obj.S(9,1);obj.S(6,1)],...
                                               "y_w",obj.water.thermalCond*obj.S(5,1)./ obj.dxThermalCond); %??????



            a1 = (obj.water.thermalCond * obj.S(5) ) ./ (obj.settings.heatTransferCoef_water * obj.S(6) * obj.dxThermalCond );
            obj.MOperatorConstPrm = ones(3,1)./(ones(3,1)+a1); 
            % a1 = obj.settings.heatTransferCoef_water * obj.S(6) * dxThermalCond + obj.water.thermalCond * obj.S(5);
            % a2 = obj.settings.heatTransferCoef_water * obj.S(6) * dxThermalCond;
            % obj.MOperatorConstPrm = a2./ a1;

            a1 = obj.settings.heatTransferCoef_water * obj.S(6) * obj.dxThermalCond;
            a2 = obj.water.thermalCond * obj.S(5);

            obj.tubeInterferrenceConstPrm = a1./ a2;

            obj.lowPassFilterTimePrm =struct("N",options.tau,"D",options.tau,"p",options.p,"p2",options.p2);

            obj.peltierPrm = obj.peltier;

        end



        
    end
end