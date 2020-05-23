%##########################################################
%第三阶段：MATLAB中进行SDWT，获得特征向量
%##########################################################
function featureFileName = SDWTgetFeature(fileName,decomLevel,K,char) 
	readfile = csvread(fileName);	
	[featureFileName,H1L1,V1L1,D1L1] = swt2(readfile,decomLevel,'haar');
	csvwrite('featureFileName.csv', featureFileName(:,(power(char,K)*2*(decomLevel-1)+1):end))
exit
end
