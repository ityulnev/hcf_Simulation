function answer=test_energy(t,f,prop_IntensityF,hcf_IntensityF,hcf_IntensityF2,prop_IntensityT,hcf_IntensityT,hcf_IntensityT2)
testCase = MyTestClass_EnergyConservation;
testCase.f=f;
testCase.t=t;
testCase.prop_IntensityF=prop_IntensityF;
testCase.prop_IntensityT=prop_IntensityT;
testCase.hcf_IntensityF=hcf_IntensityF;
testCase.hcf_IntensityT=hcf_IntensityT;
testCase.hcf_IntensityF2=hcf_IntensityF2;
testCase.hcf_IntensityT2=hcf_IntensityT2;
answer = run(testCase);
end