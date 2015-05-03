function code=generate_morse(file)
% First step: Generate 300x20 random strings and save to text file
% You can use 
% http://www.random.org/strings/?num=300&len=20&digits=on&upperalpha=on&loweralpha=off&unique=on&format=plain&rnd=new
% or
% http://www.unit-conversion.info/texttools/random-string-generator/


fid = fopen(file,"r");
fid2 = fopen("kaggle.csv","w")
 
Fs = 8000;  % 8 KHz sampling rate 
Tune = 600; % 600 Hz signal 

row = 1;
fprintf(fid2, "ID,Prediction, Usage, file,SNR,Tune,WPM\n");
while (!feof(fid))
	txt = fgetl (fid);
	fname = sprintf("cw%03d.wav",row);
	SNR = randi([-12,20]);		% SNR between -12 dB ... + 20 dB
	WPM = randi([12,60]);    	% Speed between 12 ... 60 WPM
	Set = randi([1,3]);			% Set 1=training, 2=validation, 3=testing
	%Tune = randi([300,1200]);	% Null beat between 300 ... 1200 Hz
	x = morse(txt,fname,SNR,Tune,Fs,WPM);
	fprintf(fid2, "%d,%s,%d,%s,%d,%d,%d\n",row,txt,Set,fname,SNR,Tune,WPM);
	printf("%s,row:%d,%s,SNR:%d,Tune:%d,WPM:%d\n",txt,row,fname,SNR,Tune,WPM);
	row = row + 1;
end
fclose(fid); 
fclose(fid2);
