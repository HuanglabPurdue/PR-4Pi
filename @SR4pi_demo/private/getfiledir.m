function filedir = getfiledir(filefolder,label,N)
files = dir([filefolder,'\*',label,'*.mat']);
filename = files(N).name;
filedir = fullfile(filefolder,filename);
