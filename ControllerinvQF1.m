classdef ControllerinvQF1 < Operator & Integrator & handle
    properties
    end

    properties (SetAccess = immutable, GetAccess=public)
    end

    properties (SetAccess = public, GetAccess = public) 
    end

    methods
        function obj = ControllerinvQF1(options)
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
            obj.stateVariable = zeros(4,1);
        end
        
        function dxdt = getDxdt(obj,input13,tStateVariable)
            dxdt = zeros(4,1);
            dxdt([1,3]) = -2*obj.prm.lowPassFilterTimePrm.p * tStateVariable([1,3]) + obj.prm.lowPassFilterTimePrm.p^2 *(-tStateVariable([2,4]) +input13);
            dxdt([2,4]) = tStateVariable([1,3]);
        end

        function operatorOutput =  calcNextCycle(obj,input13)
            obj.stateVariable = cMoritaRungeKuttaMethod(obj,input13);
            operatorOutput = obj.prm.lowPassFilterTimePrm.N * obj.stateVariable([1,3]) + obj.stateVariable([2,4]); 
            % operatorOutput = obj.stateVariable([2,4]);
        end

        function prevVariable = getPrevVariable(obj)
            prevVariable = obj.stateVariable;
        end

    end
end