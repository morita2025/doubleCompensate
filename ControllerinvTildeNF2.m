classdef ControllerinvTildeNF2 < Operator & Integrator & handle
    properties
    end

    properties (SetAccess = immutable, GetAccess=public)
    end

    properties (SetAccess = public, GetAccess = public) 
    end

    methods
        function obj = ControllerinvTildeNF2(options)
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
            obj.stateVariable = zeros(6,1);
        end
        
        function dxdt = getDxdt(obj,input13,tStateVariable)
            dxdt = zeros(6,1);
            dxdt([1,4]) = -3*obj.prm.lowPassFilterTimePrm.p2 * tStateVariable([1,4]) - 3*obj.prm.lowPassFilterTimePrm.p2^2 *tStateVariable([2,5])... 
            +obj.prm.lowPassFilterTimePrm.p2^3 *(-tStateVariable([3,6]) +input13);
            dxdt([2,3,5,6]) = tStateVariable([1,2,4,5]);
        end

        function operatorOutput =  calcNextCycle(obj,input13)
            obj.stateVariable = cMoritaRungeKuttaMethod(obj,input13);
            p_3 = obj.prm.A_w([1,3]);
            p_4 = obj.prm.interferrenceConstPrm.y_aw([1,3]) ./ [obj.prm.mw_cw13; obj.prm.mw_cw13];
            operatorOutput = 1./p_4 .* ( obj.prm.lowPassFilterTimePrm.D.*obj.stateVariable([1,4]) + (obj.prm.lowPassFilterTimePrm.D.*p_3 + 1).*obj.stateVariable([2,5]) + p_3.* obj.stateVariable([3,6]));
            %  = obj.prm.lowPassFilterTimePrm.N * obj.stateVariable([1,3]) + obj.stateVariable([2,4]); 
            
            % operatorOutput = obj.stateVariable([3,6]);
        end

        function prevVariable = getPrevVariable(obj)
            prevVariable = obj.stateVariable;
        end

    end
end