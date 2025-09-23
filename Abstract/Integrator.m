classdef (Abstract) Integrator  < handle
    properties
        stateVariable
        dxdtPreserved
    end

    methods
        getDxdt(obj)
    end
end