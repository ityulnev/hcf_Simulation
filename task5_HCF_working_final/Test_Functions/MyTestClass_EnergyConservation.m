classdef MyTestClass_EnergyConservation < matlab.unittest.TestCase
   properties
        t,f,prop_IntensityF,hcf_IntensityF,hcf_IntensityF2,prop_IntensityT,hcf_IntensityT,hcf_IntensityT2
        rel_tolerance=0.01;
   end

    methods (Test)
        %%Tests for Energy Conservation
        %Checks energy conservation between propagated pulse and output data both in frequency
        function test_EnergyConservation1(testCase)
            actSolution = get_Fluence(testCase.f,testCase.hcf_IntensityF2);
            expSolution = get_Fluence(testCase.f,testCase.prop_IntensityF);
            testCase.verifyEqual(round(actSolution,1),round(expSolution,1),'RelTol',testCase.rel_tolerance);
        end
        %Checks energy conservation between propagated pulse and output data both in time        
        function test_EnergyConservation2(testCase)
            actSolution = get_Fluence(testCase.t,testCase.hcf_IntensityT2);
            expSolution = get_Fluence(testCase.t,testCase.prop_IntensityT);
            testCase.verifyEqual(round(actSolution,1),round(expSolution,1),'RelTol',testCase.rel_tolerance);
        end
        
        %Checks energy conservation between input pulse and propagated pulse both in frequency
        function test_EnergyConservation3(testCase)
            actSolution = get_Fluence(testCase.f,testCase.hcf_IntensityF);
            expSolution = get_Fluence(testCase.f,testCase.prop_IntensityF);
            testCase.verifyEqual(round(actSolution,1),2.2*round(expSolution,1),'RelTol',testCase.rel_tolerance);
        end
        
        %Checks energy conservation between input pulse and propagated pulse both in time
        function test_EnergyConservation4(testCase)
            actSolution = get_Fluence(testCase.t,testCase.hcf_IntensityT);
            expSolution = get_Fluence(testCase.t,testCase.prop_IntensityT);
            testCase.verifyEqual(round(actSolution,1),2.2*round(expSolution,1),'RelTol',testCase.rel_tolerance);
        end
        
        %Checks energy conservation between input pulse in frequency and in
        %time, so the Fourier Transform
        function test_EnergyConservation5(testCase)
            actSolution = get_Fluence(testCase.f,testCase.hcf_IntensityF);
            expSolution = get_Fluence(testCase.t,testCase.hcf_IntensityT);
            testCase.verifyEqual(round(actSolution,1),round(expSolution,1),'RelTol',testCase.rel_tolerance);
        end
        
        %Checks energy conservation between output pulse in frequency and in
        %time, so the Fourier Transform
        function test_EnergyConservation6(testCase)
            actSolution = get_Fluence(testCase.f,testCase.hcf_IntensityF2);
            expSolution = get_Fluence(testCase.t,testCase.hcf_IntensityT2);
            testCase.verifyEqual(round(actSolution,1),round(expSolution,1),'RelTol',testCase.rel_tolerance);
        end
        
    end
    
end 


