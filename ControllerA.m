classdef ControllerA < Operator & Integrator & handle
    properties
    end

    properties (SetAccess = immutable, GetAccess=public)
    end

    properties (SetAccess = public, GetAccess = public) 
    end

    methods
        function obj = ControllerA(options)
            arguments                
                options.prm = [];
                options.dt  =1;
                options.cycleNum  = 1;
                options.comment = "";
            end
            obj.prm = options.prm;
            obj.dt = options.dt;
            obj.cycleNum = options.cycleNum;
            obj.comment = options.comment;

            obj.debugVectorCell = {};
            obj.stateVariable = zeros(2,1);
        end
        
        function dxdt = getDxdt(obj,input13,tStateVariable)
            dxdt = 1 ./ obj.prm.lowPassFilterTimePrm.p_A .* (-tStateVariable +input13);
        end

        function operatorOutput =  calcNextCycle(obj,ya13)
            obj.stateVariable = cMoritaRungeKuttaMethod(obj,ya13);
            operatorOutput = obj.stateVariable;
        end

    end
end