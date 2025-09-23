classdef (Abstract) Operator < handle
    properties
        prm
        dt
        comment
        cycleNum
        debugVectorCell
    end

    methods
        calcNextCycle(obj)
        getPrevVariable(obj)
    end
end