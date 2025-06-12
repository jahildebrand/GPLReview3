function sp = draw_spectogram3(data,fftl,skip,plot_length,...
    sample_freq,low_freq,high_freq,brightness,yes_Bird_labels, yes_A_labels, ...
    markers, j_included,bt,whiten,dim_coords,A_labels)
warning off;

nrec=(length(data)/fftl)*(fftl/skip); %allows for last call to be seen fully
%nrec=(sample_freq*plot_length/fftl)*(fftl/skip);

freq=[0:fftl/2]/fftl*sample_freq;
[~,low]=min(abs(freq-low_freq));
[~,high]=min(abs(freq-high_freq));
win=hamming(fftl);

%x=data;
x=[data;zeros(fftl,1)];[x1,x2]=size(x);
if x2>x1
    x=x';end

nrec = round(nrec);
sp=zeros(high-low+1,nrec);
for j=1:nrec
    start=(j-1)*skip+1;
    finish=start+fftl-1;
    q=fft(x(start:finish).*win);
    sp(:,j)=abs(q(low:high));end


if(whiten==1)

    nfre=high-low+1;
    [spc,~,mu]=whiten3(sp);
    spc=abs((spc./(mu*ones(1,nrec))).')';

    [spc,~,mu]=whiten3(spc');
    spc=abs(spc'./(ones(nfre,1)*mu'));


    cross=ones(3,3);cross(2,2)=4;cross(1,3)=0;
    cross(3,1)=0;cross(1,1)=0;cross(3,3)=0;cross=cross/8;
    spc=conv2(spc,cross,'same');
    sp=spc;
end

sz1=size(dim_coords);

if(sz1(2)>1)
    dim=zeros(size(sp));
    dim(dim_coords(1,1):dim_coords(1,2),:)=1;

    imagesc((dim.*sp).^((brightness)*2/3));axis xy;

else
    imagesc((sp).^((brightness)*2/3));axis xy;
end

[sz1,~]=size(sp);
% set(gca,'XTick',[1:floor(nrec/20):nrec])
% set(gca,'XTickLabel',[0:plot_length/20:plot_length], 'TickDir', 'out')
xticks = [1:floor(nrec/20):nrec];
xticklabels = linspace(0, plot_length, length(xticks));
set(gca, 'XTick', xticks, ...
    'XTickLabel', sprintfc('%.1f', xticklabels), ...
    'TickDir', 'out');

set(gca,'YTick',[1:20:sz1])
set(gca,'YTickLabel',[low_freq:round(20*(high_freq-low_freq)/sz1):high_freq])
%

if (yes_Bird_labels == 1)

    hold on;

    % Detect the current y-axis limits of the spectrogram
    y_limits = ylim;
    y_bottom = y_limits(1);  % Lowest visible frequency
    y_top = y_limits(2);  % Highest visible frequency

    % Compute the dynamic height as 1.5% of the visible y-range
    box_height = 0.015 * (y_top - y_bottom);

    i = 1;
    while (markers(i) < (plot_length * sample_freq) &&  j_included(i)  <= height(bt))

        % Plot vertical marker lines
        plot([(markers(i) - (fftl / 2 - skip)) / skip, ((markers(i) - (fftl / 2 - skip)) / skip)], ...
            [y_bottom, y_top], 'w');
        plot([1+((markers(i) - (fftl / 2 - skip)) / skip), 1+((markers(i) - (fftl / 2 - skip)) / skip)], ...
            [y_bottom, y_top], 'k');
        % Assign color based on boxnumber (bt column 3)
        boxnumber = bt.b_cross_flag(j_included(i)); % Use field name

        % Set color based on boxnumber
        if boxnumber == 0
            boxcolor = 'r'; % Red for 0 false detection
        elseif boxnumber == 1
            boxcolor = 'g'; % Green for 1 true detection
        else
            boxcolor = 'b'; % Default to blue for other values (optional)
        end

        % Set rectangle position exactly at the lower edge of the plot
        rectangle('Position', [(markers(i) - (fftl / 2 - skip)) / skip, ...  % X position
            y_bottom, ...  % Y position (now dynamically detected)
            ((markers(i + 1) - (fftl / 2 - skip)) / skip) - ((markers(i) - (fftl / 2 - skip)) / skip), ... % Width
            box_height], ... % Height (3% of current y-range)
            'FaceColor', boxcolor, ...
            'EdgeColor', 'k', ... % Add black border
            'LineWidth', 1.5); % Set border thickness
        i = i + 1;
        if i > length(j_included)
            break
        end
    end

    hold off;

end

if (yes_A_labels == 1)

    hold on;

    % Detect the current y-axis limits of the spectrogram
    y_limits = ylim;
    y_bottom = y_limits(1);  % Lowest visible frequency
    y_top = y_limits(2);  % Highest visible frequency

    % Compute the dynamic height as 1.5% of the visible y-range
    box_height = 0.015 * (y_top - y_bottom);

    i = 1;
    while (markers(i) < (plot_length * sample_freq) &&  j_included(i)  <= height(bt))

        % Assign color based on boxnumber (bt column 3)
        boxnumber = A_labels(j_included(i)); % Use field name

        % Set color based on boxnumber
        if boxnumber == 0
            boxcolor = 'r'; % red for no plane
        elseif boxnumber == 1
            boxcolor = 'g'; % green for jet
        elseif boxnumber == 2
            boxcolor = 'c'; % cyan for prop
        else
            boxcolor = 'b'; % Default to blue for other values (optional)
        end

        % Set rectangle position exactly at the lower edge of the plot
        rectangle('Position', [(markers(i) - (fftl / 2 - skip)) / skip, ...  % X position
            y_top - box_height, ...  % Y position (now dynamically detected)
            ((markers(i + 1) - (fftl / 2 - skip)) / skip) - ((markers(i) - (fftl / 2 - skip)) / skip), ... % Width
            y_top], ... % Height (3% of current y-range)
            'FaceColor', boxcolor, ...
            'EdgeColor', 'k', ... % Add black border
            'LineWidth', 1.5); % Set border thickness
        i = i + 1;
        if i > length(j_included)
            break
        end
    end

    hold off;

end

end
