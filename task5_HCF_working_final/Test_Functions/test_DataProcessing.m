function answer=test_DataProcessing(x_data,y_data,x_data2,y_data2,y_data2c,x_data3,y_data3,x_data3c,x_data4,y_data4)
testCase2 = MyTestClass_DataProcessing;
testCase2.x_data=x_data;
testCase2.y_data=y_data;
testCase2.x_data2=x_data2;
testCase2.y_data2=y_data2;
answer = run(testCase2);
end