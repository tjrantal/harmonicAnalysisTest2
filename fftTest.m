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

%Calculate harmonic analysis
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

%Reconstuct signal from harmonics
function reconstructed = reconstructFromAmps(harmonicAmps,period,dt)
	t = 0:dt:period;
	reconstructed = zeros(size(t,1),size(t,2));
	for i =1:length(harmonicAmps)
		reconstructed = reconstructed+harmonicAmps(i)*sin(2*pi*(1/period)*i*t);
	end
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
dt = 1/250;
stepDurations = [0.5 0.6];
stepFreqs = (1./stepDurations)/2;
fs = [1*stepFreqs; 2*stepFreqs; 3*stepFreqs; 4*stepFreqs;5*stepFreqs; 6*stepFreqs; 7*stepFreqs;8*stepFreqs; 9*stepFreqs; 10*stepFreqs;11*stepFreqs; 12*stepFreqs;];
amps = [0.05,0.03; ...	%1
		1,0.9; ...		%2
		0.11,0.12; ...	%3
		0.3,0.25; ... 	%4
		0.11,0.12; ...	%5
		0.25,0.23; ... 	%6
		0.13,0.125; ...	%7
		0.27,0.28; ... 	%8
		0.13,0.125; ...	%9
		0.23,0.22; ... 	%10
		0.08,0.085; ...	%11
		0.16,0.16; ... 	%12
		];
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
%gTK = 'fltk';
graphics_toolkit(gTK);

for sit = 1:length(sError)
	if showFig
		titles = {sprintf("signal serr %.2f",sError(sit)),'hanSignal','ffAmp','hanFftAmp','harmonicAnalysis','hanHarmonicAnalysi','reco','hanReco'};
		fh = figure('__graphics_toolkit__',gTK,'position',[10 10 1000 1000]);
		for i = 1:8
			ax(i) = subplot(4,2,i);
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
		reconstructed = reconstructFromAmps(harmonicAmps,strideDuration,dt);
		hHarmonicAmps = harmonicAmp(hanSig,t,strideDuration+sError(sit));
		hReconstructed = reconstructFromAmps(2*hHarmonicAmps,strideDuration,dt);
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
			set(gca,'xlim',[0 3.5],'ylim',[-1.5 2])
			set(fh,'currentaxes',ax(2));
			plot(t,hanSig,'color',cMap(it,:));
			set(gca,'xlim',[0 3.5],'ylim',[-1.5 2])
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
			set(fh,'currentaxes',ax(7));
			plot(t(1:length(reconstructed)),reconstructed,'color',cMap(it,:));
			set(gca,'xlim',[0 3.5],'ylim',[-1.2 2])
			set(fh,'currentaxes',ax(8));
			plot(t(1:length(hReconstructed)),hReconstructed,'color',cMap(it,:));
			set(gca,'xlim',[0 3.5],'ylim',[-1.2 2])
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
