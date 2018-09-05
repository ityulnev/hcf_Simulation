classdef SolverTest < matlab.unittest.TestCase

    
    methods (Test)
        function testRealSolution(testCase)
            actSolution = doubler(8);
            expSolution = 16;
            testCase.verifyEqual(actSolution,expSolution);
        end
    end
    
end 


function out = doubler(in)
out = 2*in;
end