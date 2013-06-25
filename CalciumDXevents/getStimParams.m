function [stimuli] = getStimParams(fhandle, channelToFilter, times, stim_dur, stim_ISI, first_stim_tim, number_of_stim)
%stim_dur %stimulus duration in seconds
%If the Stimulus was marked with an Event marker as a channel in the file then use the following:
% [stimulusParams] = getStimParams(2,3,times,0.1)
%If the Precise times of the Events are known then leave the channelToFilter empty, and instead include the interstimulus interval (stim_ISI in seconds), time of first stimulus (first_stim_tim), and number of stimuli (number_of_stim), like if the stimulus waveform software was triggered off the start of the recording (Prairie View Imaging) software.
%[stimulusParams] = getStimParams(2,[],times,0.1,30,10,10) 
%James B. Ackman (c) 7/27/2011 

[fhandle, channels]=scParam(fhandle);
% data=(channels{channelToFilter}.adc(:,1))';
% stim_dur=0.3;  %LED stim duration in seconds
stim_dur=stim_dur*1e06; %convert to microseconds

if ~isempty(channelToFilter)
    for numStim = 1:numel(channelToFilter)
        chan = channelToFilter(numStim);
        for i=1:numel(channels{chan}.tim(:))
            % channelToFilter=3;
            % i=11;
            
            if stim_dur < min(diff(times))
                idx=find(times < (channels{chan}.tim(i)+stim_dur));
                frame_indices = idx(end);
            else
            frame_indices=find(channels{chan}.tim(i) < times & times < (channels{chan}.tim(i)+stim_dur));
            end
            
            if ~isempty(frame_indices)
                %             if times(frame_indices(1)-1) < channels{channelToFilter}.tim(i) && channels{channelToFilter}.tim(i) < times(frame_indices(1))-1.235 %minus the small interframe interval which varies around 1.235msec to 1.60msec, even though in this line it's wrong anyways (should be conv to microseconds)
                if times(frame_indices(1)-1) < channels{chan}.tim(i) && channels{chan}.tim(i) < times(frame_indices(1))
                    frame_indices=[frame_indices(1)-1; frame_indices];
                end
                stimuli{numStim}.stimulusParams{i}.frame_indices=frame_indices;
                stimuli{numStim}.stimulusParams{i}.frame_times=times(frame_indices);
                stimuli{numStim}.stimulusParams{i}.stimulus_times=[channels{chan}.tim(i) channels{chan}.tim(i)+stim_dur];
            end
        end
        %uidialog for stim description
        prompt = {'description'};
        dlg_title = ['Stimulus desc for ' 'ch' num2str(chan) ', 1st stim at ' num2str(stimuli{numStim}.stimulusParams{i}.stimulus_times(1) * 1e-06) 's'];
        num_lines = 1;
        def = {'LED.stim.R.eye'};
        answer = inputdlg(prompt,dlg_title,num_lines,def);
        stimuli{numStim}.description = answer;
    end
else
    tim=linspace(first_stim_tim,(number_of_stim-1)*stim_ISI+first_stim_tim,number_of_stim);
    tim=tim.*1e06;
    
    for i=1:numel(tim)
        frame_indices=find(tim(i) < times & times < (tim(i)+stim_dur));
        if ~isempty(frame_indices)
            %if times(frame_indices(1)-1) < tim(i) && tim(i) < times(frame_indices(1))-1.235  %minus the small interframe interval which varies around 1.235msec to 1.60msec, even though in this line it's wrong anyways (should be conv to microseconds)
            if times(frame_indices(1)-1) < tim(i) && tim(i) < times(frame_indices(1))
                frame_indices=[frame_indices(1)-1; frame_indices];
            end
            stimuli{1}.stimulusParams{i}.frame_indices=frame_indices;
            stimuli{1}.stimulusParams{i}.frame_times=times(frame_indices);
            stimuli{1}.stimulusParams{i}.stimulus_times=[tim(i) tim(i)+stim_dur];
        end
    end
    %uidialog for stim description
    prompt = {'description'};
    dlg_title = 'Stimulus description';
    num_lines = 1;
    def = {'LED.stim.R.eye'};
    answer = inputdlg(prompt,dlg_title,num_lines,def);
    stimuli(1).description = answer;
end