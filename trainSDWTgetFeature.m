function featureFileName = trainSDWTgetFeature(fileName,decomLevel,K) 
	readfile = csvread(fileName);	
	[featureFileName,H1L1,V1L1,D1L1] = swt2(readfile,decomLevel,'haar');
	csvwrite('trainfeatureFileName.csv', featureFileName(:,(power(4,K)*2*(decomLevel-1)+1):end));
exit
end
