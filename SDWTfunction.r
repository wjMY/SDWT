#########################################################
#Sequence reading and stitching
#######################################################
DNAalphabet="ACGT"
proteinalphabet="ACDEFGHIKLMNPQRSTVWY"

Connect=function(x,n1)
{
	y=array()
	m=0
	i=1
	for(j in 1:n1)
	{
		temp_x=NULL
		if(substr(x[i],1,1)==">")
		{#If the first character read is ">", i adds 1 and reads the next line
			m=m+1
			i=i+1
		}	
		while(m==j)
		{#When m=j, control the spliced ​​newline
			if(substr(x[i],1,1)!=">"&i!=(length(x)+1))
			{#Splicing if the read first character is not ">"
				temp_x=paste(temp_x,x[i],sep="")
				i=i+1
			}
			else break#Otherwise jump out of the loop without stitching, start the next stitching
		}
		y[j]=temp_x	
	}
	return(y)
}
#######################################################
#Read the HOG file to complete the stitching
#################################################
readFile = function(filename)
{
	read_file = read.table(filename,header = FALSE)
	family_nums = nrow(read_file)
	seqs = NULL
	seqLabel = NULL
	label = 0
	path = strsplit(filename,split="/")[[1]][1]
	for(i in 1:family_nums)
	{
		temp_a=readLines(paste(path,"/",read_file[i,1],sep=""),encoding = "UTF-8")#读取文件中每个family的序列
		n=0
		label = label+1
		for(j in 1:length(temp_a))
		{
			if(substr(temp_a[j],1,1)==">")
			n=n+1
		}
		temp_label = rep(label,n)
		seqLabel = c(seqLabel,temp_label)
		oneSeq_connet = NULL
		oneSeq_connet = Connect(temp_a,n)#Splicing the sequence in the current family	
		seqs=c(seqs,oneSeq_connet)	#Combine the sequence from each family read in one file
	}
	outputSeqLabel = list(seqs=seqs,seqLabel=seqLabel)
	return (outputSeqLabel)
}
#######################################################
#Read the HOG file to complete the stitching
#################################################
testReadFile = function(filename)
{
	read_file = read.table(filename,header = FALSE)
	family_nums = nrow(read_file)
	allHOGseqs = NULL
	path = strsplit(filename,split="/")[[1]][1]
	for(i in 1:family_nums)
	{
		temp_a=readLines(paste(path,"/",read_file[i,1],sep=""),encoding = "UTF-8")#读取文件中每个family的序列
		n=0
		for(j in 1:length(temp_a))
		{
			if(substr(temp_a[j],1,1)==">")
			n=n+1
		}
		oneSeq_connet = NULL
		oneSeq_connet = Connect(temp_a,n)#Splicing the sequence in the current family
		allHOGseqs=c(allHOGseqs,oneSeq_connet)	#Combine the sequence from each family read in one file
	}
	return (allHOGseqs)
}

#2#############################################################
#Statistical word frequency, can only operate in Linux environment
kmer<-function(seq, k, alphabet, len=length(seq))
{
    dyn.load("Rkmer.so")
	alphasize=nchar(alphabet)
	matrix(.C("Rkmer",seq,as.integer(k),as.integer(len),alphabet,retArray=integer(len*alphasize^k))$retArray,nrow=len,ncol=alphasize^k,byrow=T)
}

TrainFeature =function(K,dat,decomLevel,char,alphabet)
{
	angle=360/char^K
	K_word=array(0:0,c(char^K,2))
	angle_sum = 0
	for(i in 1:char^K)
	{
		angle_sum = angle_sum + angle
		K_word[i,1]=sinpi(angle_sum /180)
		K_word[i,2]=cospi(angle_sum/180)	
	}
	##########
	wordFrequence = kmer(dat,K,alphabet,length(dat))
	StandFrequency = scale(wordFrequence)
	means = attr(StandFrequency ,"scaled:center")
	SD = attr(StandFrequency ,"scaled:scale")
	fre_real=array(0:0,c(nrow(StandFrequency),char^K))
	fre_imag=array(0:0,c(nrow(StandFrequency),char^K))
	for (i in 1:nrow(StandFrequency)) 
	{  
		#Word Frequency Multiplied by Complex Numbers
		for(j in 1:char^K)
		{
			fre_real[i,j]=StandFrequency[i,j]*(K_word[j,1])#Real
			fre_imag[i,j]=StandFrequency[i,j]*(K_word[j,2])#Imaginary
		}
	}
	complexFileName = cbind(fre_real,fre_imag)#Real and imaginary mergers
	if (nrow(complexFileName)%%(2^decomLevel)!=0)
	{
		add = array(0:0,c((2^decomLevel)-nrow(complexFileName)%%(2^decomLevel),ncol(complexFileName)))
		complexFileName = rbind(complexFileName,add)
	}
	write.table(complexFileName,file="traincomplexFileName.csv",append=TRUE,row.names = FALSE, col.names=FALSE, sep=",")
	mat=paste("matlab -nojvm -nodisplay -nosplash -nodesktop -r \"trainSDWTgetFeature(\'traincomplexFileName.csv\'",decomLevel,K,sep=",")
	mat=paste(mat,")\"",sep="")
	system(mat)
	trainseqfeature = read.csv("trainfeatureFileName.csv",header=FALSE)[1:length(dat),]
	output = list(trainseqfeature=trainseqfeature,means=means,SD=SD)
	return (output)
}
#######################################################
TestFeature = function(K,dat,decomLevel,center,sde,char,alphabet)
{
	angle=360/char^K
	K_word=array(0:0,c(char^K,2))
	angle_sum = 0
	for(i in 1:char^K)
	{
		angle_sum = angle_sum + angle
		K_word[i,1]=sinpi(angle_sum /180)
		K_word[i,2]=cospi(angle_sum/180)	
	}
	##########
	wordFrequence = kmer(dat,K,alphabet,length(dat))
	StandFrequency=array(0:0,c(nrow(wordFrequence),char^K))
	for (i in 1:nrow(wordFrequence))
	{
		for(j in 1:char^K)
		{
			StandFrequency[i,j]=(wordFrequence[i,j]-center[j])/sde[j]
		}
	}
	fre_real=array(0:0,c(nrow(StandFrequency),char^K))
	fre_imag=array(0:0,c(nrow(StandFrequency),char^K))
	for (i in 1:nrow(StandFrequency)) 
	{  
		for(j in 1:char^K)
		{
			fre_real[i,j]=StandFrequency[i,j]*(K_word[j,1])#
			fre_imag[i,j]=StandFrequency[i,j]*(K_word[j,2])#
		}
	}
	complexFileName = cbind(fre_real,fre_imag)#
	if (nrow(complexFileName)%%(2^decomLevel)!=0)
	{
		add = array(0:0,c((2^decomLevel)-nrow(complexFileName)%%(2^decomLevel),ncol(complexFileName)))
		complexFileName = rbind(complexFileName,add)
	}
	write.table(complexFileName,file="testcomplexFileName.csv",append=TRUE,row.names = FALSE, col.names=FALSE, sep=",")
	mat=paste("matlab -nojvm -nodisplay -nosplash -nodesktop -r \"testSDWTgetFeature(\'testcomplexFileName.csv\'",decomLevel,K,sep=",")
	mat=paste(mat,")\"",sep="")
	system(mat)
	testfeature = read.csv("testfeatureFileName.csv",header=FALSE)[1:length(dat),]
	return (testfeature)
}

Feature = function(K,dat,decomLevel,char,alphabet)
{
	angle=360/char^K
	K_word=array(0:0,c(char^K,2))
	angle_sum = 0
	for(i in 1:char^K)
	{
		angle_sum = angle_sum + angle
		K_word[i,1]=sinpi(angle_sum /180)
		K_word[i,2]=cospi(angle_sum/180)	
	}
	##########
	wordFrequence = kmer(dat,K,alphabet,length(dat))
	StandFrequency = scale(wordFrequence)
	
	fre_real=array(0:0,c(nrow(StandFrequency),char^K))
	fre_imag=array(0:0,c(nrow(StandFrequency),char^K))
	for (i in 1:nrow(StandFrequency)) 
	{ 
		for(j in 1:char^K)
		{
			fre_real[i,j]=StandFrequency[i,j]*(K_word[j,1])#
			fre_imag[i,j]=StandFrequency[i,j]*(K_word[j,2])#
		}
	}
	complexFileName = cbind(fre_real,fre_imag)#
	if (nrow(complexFileName)%%(2^decomLevel)!=0)
	{
		add = array(0:0,c((2^decomLevel)-nrow(complexFileName)%%(2^decomLevel),ncol(complexFileName)))
		complexFileName = rbind(complexFileName,add)
	}
	write.table(complexFileName,file="complexFileName.csv",append=TRUE,row.names = FALSE, col.names=FALSE, sep=",")
	mat=paste("matlab -nojvm -nodisplay -nosplash -nodesktop -r \"SDWTgetFeature(\'complexFileName.csv\'",decomLevel,K,char,sep=",")
	mat=paste(mat,")\"",sep="")
	system(mat)
	feature = read.csv("featureFileName.csv",header=FALSE)[1:length(dat),]
	return (feature)
}


