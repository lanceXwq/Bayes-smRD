classdef ExperimentalData < handle

    properties
        units
        num_roi uint32 = 1
        param = ExperimentalParam()
        detn = ExperimentalDetn()
    end

    methods

        function obj = ExperimentalData(varargin)

            if nargin == 1
                num_roi = varargin{1};
                obj.detn(num_roi, 1) = ExperimentalDetn();
                obj.num_roi = num_roi;
            end

        end

    end

end
