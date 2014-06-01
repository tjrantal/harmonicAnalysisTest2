%Create two steps, one longer, and one shorter to create a full
%stride. The stride is then repeated 2.5 to 3 times

close all;
clear all;
clc;

function fftAmps = fftAmp(signal) 
	fftLength = double(int32(length(signal)/2));
	fftSignal = fft(signal);
	fftSignal= fftSignal./(fftLength+1);
    fftSignal([1 fftLength+1]) = fftSignal([1 fftLength+1])/2;
 	fftAmps = abs(fftSignal(1:(fftLength+1))); %Ignore the 	second half of the fft
endfunction;

function harmonicAmps = harmonicAmp(signal,t,period)
	A = ones(length(signal),1);
	if size(t,1) < size(t,2)
		t=t';
	end
	for i = 1:20
		A = [A, cos(i*2*pi*t/period),sin(i*2*pi*t/period)]; 
	end
	if size(signal,1)<size(signal,2)
		signal = signal';
	end
	b = A\signal;
	coeffs(:,1) = b([2:2:40]);
	coeffs(:,2) = b([3:2:41]);
	harmonicAmps = sqrt(sum(coeffs.^2,2));
endfunction;

%Open results file, and write a header
testResults = fopen('fftResulFile.xls',"w");
		fprintf(testResults,"sErr\toffs\tHR\thanHR\t");
		for ha = 1:20
			fprintf(testResults,"HarmonicAmp %d\t",ha);
		end
		for ha = 1:20
			fprintf(testResults,"HanHarmonicAmp %d\t",ha);
		end
		fprintf(testResults,"\n");

%Create the stride
dt = 0.005;
stepDurations = [0.6 0.5];
stepFreqs = (1./stepDurations)/2;
fs = [stepFreqs; 2*stepFreqs; 3*stepFreqs; 4*stepFreqs;];
amps = [0.1,0.08; 1,1.02; 0.11,0.12; 0.3,0.25];
stride = [];
%concat steps to a stride
for s = 1:length(stepDurations)
	t = 0:dt:(stepDurations(s)-dt);
	step = zeros(size(t));
	for i = 1:size(fs,1)
		step = step+sin(2*pi*fs(i,s)*t)*amps(i,s);
	end
	stride = [stride step];
end
timeInstances = ([1:length(stride)]-1)*dt;
%figure
%plot(timeInstances,stride)

%Create 3 full strides
origSignal = [];
for i =1:3
	origSignal = [origSignal stride];
end
%figure
%plot(origSignal)
%keyboard;	
strideDuration = sum(stepDurations);
strideFreq = 1/strideDuration;
sError = [-0.15:0.05:0.15];
sampleStop = strideDuration*1;

iterations = 10;
offset = strideDuration/2/iterations;
cMap = jet(iterations);
showFig = 1;
gTK = 'gnuplot';
graphics_toolkit(gTK);

for sit = 1:length(sError)
	if showFig
		titles = {sprintf("signal serr %.2f",sError(sit)),'hanSignal','ffAmp','hanFftAmp','harmonicAnalysis','hanHarmonicAnalysi'};
		fh = figure('__graphics_toolkit__',gTK,'position',[10 10 1000 500]);
		for i = 1:6
			ax(i) = subplot(3,2,i);
			hold on;
			title(titles{i});
		end
	end
	for it = 1:iterations
		
		t = 0:dt:(3*sampleStop-((it-1)*offset))-dt;
		%t = 0:dt:sampleStop;
		signal = origSignal(1:length(t));
		%FFT analysis
		fftLength = double(int32(length(signal)/2));
		freq= [0:fftLength]/dt/2/fftLength;
		hf = find(freq >= 20,1,'first');
		amp = fftAmp(signal);
		%Hann windowed
		hanSig = signal.*hanning(length(signal))';
		hAmp = fftAmp(hanSig);
	
		%Harmonic analysis
		harmonicAmps = harmonicAmp(signal,t,strideDuration+sError(sit));
		hHarmonicAmps = harmonicAmp(hanSig,t,strideDuration+sError(sit));
		%calc harmonic ratios
		hr = sum(harmonicAmps(2:2:20))/sum(harmonicAmps(1:2:19));
		hhr = sum(hHarmonicAmps(2:2:20))/sum(hHarmonicAmps(1:2:19));
		fprintf(testResults,"%.2f\t%.3f\t",sError(sit),t(1));
		fprintf(testResults,"%f\t%f\t",hr,hhr);	
		for ha = 1:length(harmonicAmps)
			fprintf(testResults,"%f\t",harmonicAmps(ha));
		end
		for ha = 1:length(harmonicAmps)
			fprintf(testResults,"%f\t",hHarmonicAmps(ha));
		end
		fprintf(testResults,"\n");
		if showFig
			set(fh,'currentaxes',ax(1));
			plot(t,signal,'color',cMap(it,:));
			set(fh,'currentaxes',ax(2));
			plot(t,hanSig,'color',cMap(it,:));
			set(fh,'currentaxes',ax(3));
			plot(freq(1:hf),amp(1:hf),'color',cMap(it,:))
			set(fh,'currentaxes',ax(4));
			plot(freq(1:hf),2*hAmp(1:hf),'color',cMap(it,:));
			set(fh,'currentaxes',ax(5));
			bar(harmonicAmps,'facecolor','none','edgecolor',cMap(it,:),'linewidth',3)
			set(gca,'ylim',[0 1])
			set(fh,'currentaxes',ax(6));
			bar(2*hHarmonicAmps,'facecolor','none','edgecolor',cMap(it,:),'linewidth',3)
			set(gca,'ylim',[0 1])
			%drawnow();
		end
	end
	if showFig
		%keyboard;
		drawnow();
		print('-dpng','-r200','-S2000,1000',sprintf("./figures/serr%.0f.png",sError(sit)*100));
		%keyboard;
		close;
	end
	%keyboard;
end
fclose(testResults);
