# SDWT
 A new sequence similarity analysis method based on stationary discrete wavelet transform 


1. Need to install class classification of R packages and matlab
2. Usage
a. Read file function: Reads sequence data according to the parameter familylist and belongs to the family
Function prototype : seqsAndLabel = readFile("familylist")
input: familylist, This parameter is the input file name, which contains the family file name. Each file contains the corresponding gene sequence.
output: seqsAndLabel, When the training set file is read in the classification, the returned output variable includes two columns of vectors, which are two types of different data. The first column is a training set sequence, and the second column is a training set sequence tag
When reading the test set file in the classification, the returned output variable is just the test set sequence.
When the cluster reads the file, the returned output variable contains two columns of vectors, which are two different types of data. The first column is the sequence set, and the second column is the number of corresponding sequences in each family.
b. Feature Constructor: Constructs SDWT features based on the input sequence
Function prototype: feature = Feature(K,dat,decomLevel,char,alphabet)
input:
    K: the parameter of k-mer
    dat: the vector of sequences
    decomLevel: decomposition level of wavelet transform
	char: the number of character classes contained in the sequence
	alphabet: characters included in the sequence
output:  
    featrue: the matrix of feature
    

c. Training Feature Constructor: Construct SDWT features and training data features based on the input sequence
Function prototype: traindata = TrainFeature(K,dat,decomLevel,char,alphabet)
input:
    K: the parameter of k-mer
    dat: the vector of train sequences
    decomLevel: decomposition level of wavelet transform
	char: the number of character classes contained in the sequence
	alphabet: characters included in the sequence
output: 
	trainfeature: the matrix of train feature
    means: average of k-mer
    SD: Standard deviation of k-mer.

d. Predictive Feature Constructors: Construct SDWT features and training data features based on the input sequence
Function prototype: testfeature = TestFeature(K,dat,decomLevel,center,sde,char,alphabet)
input:
    K: the parameter of k-mer
    dat: the vector of train sequences
    decomLevel: decomposition level of wavelet transform
	char: the number of character classes contained in the sequence
	alphabet: characters included in the sequence
output:  
    feature: the matrix of test feature

Usage:
Classification example:
> source("SDWTfunction.R")
> seqsAndLabel = readFile("traindata/file.lst")
> testseqs = testReadFile("testdata/file.lst")
> traindata = TrainFeature(3,seqsAndLabel$seqs,3,4,DNAalphabet)
> testfeature = TestFeature(3,testseqs,3,traindata$means,traindata$SD,4,DNAalphabet)
> library(class)
> predictResult = knn(traindata$trainseqfeature,testfeature,seqsAndLabel$seqLabel,k=1)


Clustering example:

> source("SDWTfunction.R")
> seqsAndLabel = readFile("cluster/file.lst")
> feature = Feature(3,seqsAndLabel$seqs,3,4,DNAalphabet)
> KmeansResult = kmeans(feature,100)


