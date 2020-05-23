function featureFileName = testSDWTgetFeature(testfileName,decomLevel,K) 
	readfile = csvread(testfileName);	
	[featureFileName,H1L1,V1L1,D1L1] = swt2(readfile,decomLevel,'haar');
	csvwrite('testfeatureFileName.csv', featureFileName(:,(power(4,K)*2*(decomLevel-1)+1):end))
	%csvwrite('testfeatureFileName.csv', featureFileName)
exit
end