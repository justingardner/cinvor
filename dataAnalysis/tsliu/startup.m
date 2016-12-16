% addpath('~/matlab/libsvm-3.1/matlab');
pathNames ={'mgl', 'mrTools', 'gru','CircStat2011f','/Users/gru/data/cinvor2'};
addpath('gru/mac_helpers');

for i = 1:length(pathNames)
	  if isdir(pathNames{i})
	  addpath(genpath_exclude(pathNames{i},{'.git','.svn'}));
  end
end


mrSetPref('verbose','No');
mrSetPref('site','MSU');
mrSetPref('niftiFileExtension','.img');
mrSetPref('interpMethod','nearest');
mrSetPref('overwritePolicy','Ask');
mrSetPref('importROIPath','fmri/anat/');
mrSetPref('magnet',{'GE signa 3T', 'Siemens 3T', 'Philips 3T', 'other'} ); 
mrSetPref('coil',{'head array', 'other'} ); 
mrSetPref('maxBlockSize', 250000000*5); 


fshome = getenv('FREESURFER_HOME');
fsmatlab = sprintf('%s/matlab',fshome);
if (exist(fsmatlab) == 7)
    path(path,fsmatlab);
end
clear fshome fsmatlab;
fsfasthome = getenv('FSFAST_HOME');
fsfasttoolbox = sprintf('%s/toolbox',fsfasthome);
if (exist(fsfasttoolbox) == 7)
    path(path,fsfasttoolbox);
end
clear fsfasthome fsfasttoolbox;
