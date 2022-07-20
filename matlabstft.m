%%Aquire data from antenna
ai = analoginput('winsound'); %se�al de entrada
addchannel(ai,1:2);
set(ai,'SampleRate',110250)
set(ai,'SamplesPerTrigger',110250)
disp('Recording...');
start(ai);
data = getdata(ai);
x = data(:); 
