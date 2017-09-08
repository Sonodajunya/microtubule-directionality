function mtDirectionality

% PURPOSE:
%   Computes relative microtubule growth direction from data previously analysed using
%   plusTipTracker software from the Danuser lab.
%   It outputs raw data as well as a number of figures and plots
%
% INPUT:  
%   no input required to run the function
%   However, the code will require a projlist.mat file and will go through all
%   projects within this file.
%
%   It is assumed that the projData.mat file is contained in the meta subfolders in
%   all projects, i.e. all projects have been processed by plustipgettracks already.
%
%   The user will be asked for each project to draw the major cell axis relative to which
%   angles of microtubule growth will be calculated. Confirm line with
%   doubleclick and wait until tracks plot into that figure. This takes
%   about 30s, but will vary depending upon available memory.
%
% OUTPUT:
%   for each project into the analysis directory specified in projlist.mat
%       (1) plots tracks with all segments colour-coded by direction relative to the main cell axis overlayed to first image,
%       (2) plots an angle histogram using 6º bin size
%       (3) exports data tables for histograms and raw angle data as text files
%       (4) plots kuiper statistics relative to random distribution.
%
%
% Copyright: Straube lab, University of Warwick, March 2015
% shared under New BSD license as specified below
%
% Code was tested under Matlab releases 2012a and 2012b


minTrackLength=3;

[FileName,PathName]=uigetfile('*.mat','Please choose projlist.mat file to load');
load([PathName FileName]);

for k=1:size(projList,1);

datafilename = strcat(projList(k,1).anDir ,'/roiYX.mat');
load(datafilename);

datafilename = strcat(projList(k,1).anDir ,'/meta/projData.mat');
load(datafilename);
        
imscale = projData.pixSizeNm/1000; % Microns per pixel

listOfImages = dir(projData.imDir); % Get the first image in the 'imDir' directory

ROI = imread([projData.imDir '/' listOfImages(3).name]);

ROI = im2uint16(ROI);

figure(1);
clf;
colormap(gray);
% Display the image, scaled to microns,
cellimg=imagesc([0 size(ROI,2)*imscale],[0 size(ROI,1)*imscale],ROI);
% imagesc is the same as imashow but automatically adjusts the contrast
% according to the highest and lowest values present.

figure(1);
hold on;
plot(imscale*roiYX(:,2)',imscale*roiYX(:,1)','r')
xlabel('X coords (\mum)')
ylabel('Y coords (\mum)')
axis equal

title('Please select the long axis of the cell');
[X,Y]=getline;
cellAxisVec=[X(2)-X(1) Y(2)-Y(1)];
plot([X(1) X(2)],[Y(1) Y(2)],'c-');

%%%%%%%%%%%%%%%%%%%%%%%%%%
% Extracts growth direction and computes relative angles
%%%%%%%%%%%%%%%%%%%%%%%%%%


angle=[];
tracks=zeros(size(projData.xCoord,1),3);
cols=size(projData.xCoord,2);
for i=1:size(projData.xCoord,1)
    % go through each track and keep increasing time-point (moving to next
    % column) until find a non-nan element (start of a track)
    startTrack=1;
    while startTrack<=cols && isnan(projData.xCoord(i,startTrack))
        startTrack=startTrack+1;
    end
    
    % if a track is all nans then skip the rest of the loop and go to the
    % next track (i+1)
    if startTrack>cols
        continue
    end

    % once the start of the track is found it moves through the following
    % columns until it finds the first next nan. If k is the position of
    % the first nan element after the start of the track then k-1 is the
    % end of the track.
    
    endTrack=startTrack;
    while endTrack+1<=cols && ~isnan(projData.xCoord(i,endTrack+1))
        endTrack=endTrack+1;
    end

    % Threshholds for tracks that are too small
    if endTrack-startTrack+1<minTrackLength
        continue
    end
   
   trackVec=[projData.xCoord(i,endTrack)-projData.xCoord(i,startTrack) projData.yCoord(i,endTrack)-projData.yCoord(i,startTrack)];
   angle(i)=acos(dot(cellAxisVec,trackVec)/(norm(cellAxisVec)*norm(trackVec)));
   angle(i)=180*(angle(i)/pi);
   tracks(i,:)=[i startTrack endTrack];
   
   % check if going left or right with respect to axis
   normal=cross([cellAxisVec 0],[trackVec 0]);
   if normal(3)<0
       angle(i)=-angle(i);
   end
   
   % set colours relative to growth direction, here we use four sectors of 90º each, centered on the main axis as defined by the user
   % colours are green for ±45º from the main axis (forward), blue between 45º and 135º (right), yellow between -45º and -135º (left)
   % and the remaining red (backwards).
   
   
   if abs(angle(i))>=135
       cm = ('r-');
   elseif abs(angle(i))<=45
       cm = ('g-');
   elseif angle(i)>45
       cm = ('b-');
   else
       cm = ('y-');
   end
   
   %plot each segment of the track with the colour that encodes direction
   %from start to end of that track
   
   for f=startTrack:endTrack-1
            x1 = projData.xCoord(i,f)* imscale;
            x2 = projData.xCoord(i,f+1)* imscale;
            y1 = projData.yCoord(i,f)*imscale;
            y2 = projData.yCoord(i,f+1)*imscale;
            plot([x1 x2],[y1 y2],cm);
   end
   
end
 
title('growing microtubule tracks and chosen ROIs');

    saveas(gcf,[projList(k,1).anDir '/tracks.fig']);
    saveas(gcf,[projList(k,1).anDir '/tracks.pdf']);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Angle histogram and data export for each project
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
figure(2);
clf;
set(gca,'FontSize',8);
bins=[-177:6:177];
n=hist(angle,bins);
n=n/length(angle);
bar(bins,n)

xlabel('Angles (degrees)');
ylabel('Fraction');
title(['angle distribution' projList(k,1).anDir]);
ylim([0 0.2]);
xlim([-180 +180]);

b=[bins' n'];

    dlmwrite([projList(k,1).anDir '/anglehist.txt'],b,'delimiter','\t');
    saveas(gcf,[projList(k,1).anDir '/angles.fig']);
    saveas(gcf,[projList(k,1).anDir '/angles.pdf']);


angle=angle';

save([projList(k,1).anDir '/angles.mat'],'angle');
dlmwrite([projList(k,1).anDir '/angles.txt'],angle,'delimiter', '\t');


%%%%%%%%%%%%%%%%%%%%
% Kuiper statistics for each project
%%%%%%%%%%%%%%%%%%%%

x = angle;
n  = numel(x);
x = reshape(x,1,n); % ensure x is a row vector
yn = (0:n)/n;
y2 = reshape(repmat(yn,2,1),1,2*n+2); % cumulative probability
x2 = reshape(repmat(sort(x),2,1),1,2*n); % sort and duplicate variates
x2 = [-Inf x2 Inf];
yc = x2/360 + 0.5;
Kstat = max(y2(2:end-1)-yc(2:end-1)) + max(yc(2:end-1)-y2(2:end-1));
    
figure(3);
clf;


% Create axes
axes1 = axes('YTick',[0 0.2 0.4 0.6 0.8 1],...
    'XTick',[-180 -90 0 90 180],...
    'TickDir','out',...
    'FontSize',20,...
    'FontName','Arial');
axis([-180 180 0 1])
hold(axes1,'all');

plot(x2,y2,'r')
plot(x2,yc,'b')
hold off
legend('data','random', 'Location','NorthWest')    
% Create textbox
annotation(gcf,'textbox',...
    [0.7 0.24 0.15 0.12],...
    'String',{'Kstat =' Kstat},'FontSize',20,...
    'FontName','Arial');

saveas(gcf,[projList(k,1).anDir '/kuiper.fig']);
saveas(gcf,[projList(k,1).anDir '/kuiper.pdf']);



end
    
    
end




% Copyright (c) 2015, Straube lab, University of Warwick
% All rights reserved.
% 
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are met:
%     * Redistributions of source code must retain the above copyright
%       notice, this list of conditions and the following disclaimer.
%     * Redistributions in binary form must reproduce the above copyright
%       notice, this list of conditions and the following disclaimer in the
%       documentation and/or other materials provided with the distribution.
%     * Neither the name of the organization nor the
%       names of its contributors may be used to endorse or promote products
%       derived from this software without specific prior written permission.
% 
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
% ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
% WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
% DISCLAIMED. IN NO EVENT SHALL <COPYRIGHT HOLDER> BE LIABLE FOR ANY
% DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
% (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
% LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
% ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
% (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
% SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.