classdef MyTestClass_DataProcessing < matlab.unittest.TestCase
   properties
        x_data,y_data,x_data2,y_data2,y_data2c,x_data3,y_data3,x_data3c,x_data4,y_data4
        rel_tolerance=0.01;
   end

    methods (Test)
        %%Tests for data processing
        function test_DataConversion(testCase2)
            actSolution = get_Fluence(testCase2.x_data,testCase2.y_data);
            expSolution = get_Fluence(testCase2.x_data2,testCase2.y_data2);
            testCase2.verifyEqual(round(actSolution,1),round(expSolution,1),'RelTol',testCase2.rel_tolerance);
        end
        
    end
    
end 