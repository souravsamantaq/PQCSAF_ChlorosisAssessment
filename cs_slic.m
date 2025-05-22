function [l,Sp,bestnest]=cs_slic(im)
%if nargin<1,
% Number of nests (or different solutions)0013_0014_LEAF
n=10;
%end





im=imresize(im,1);

[h,w,d]=size(im);
Lfitness=zeros(n,h,w);
gbestL=zeros(h,w);

% Discovery rate of alien eggs/solutions
pa=0.25;

%% Change this if you want to get better results
% Tolerance
%Tol=1.0e-5;
Total_itr=10;
%% Simple bounds of the search domain
% Lower bounds
nd=3; 
% Lb=-5*ones(1,nd); 
% % Upper bounds
% Ub=5*ones(1,nd);
%--------------range-------------
Lb(1,1)=300;Lb(1,2)=10;Lb(1,3)=1;
Ub(1,1)=500;Ub(1,2)=40;Ub(1,3)=4;
% Lb(1,1)=50;Lb(1,2)=30;Lb(1,3)=2;
% Ub(1,1)=100;Ub(1,2)=40;Ub(1,3)=4;


% Random initial solutions
for i=1:n
    for k=1:nd
     nest(i,k)=Lb(k)+(Ub(k)-Lb(k))*rand();
    end
end

% Get the current best
%fitness=10^10*ones(n,1);
fitness=50000*ones(n,1);
[fmin,bestnest,nest,fitness,Lbest]=get_best_nest(im,nest,nest,fitness,Lfitness);
gbestL=Lbest;
N_iter=1;
%% Starting iterations
tic
while (N_iter<=Total_itr)%(fmin>Tol),

    % Generate new solutions (but keep the current best)
     [nest]=get_cuckoos(nest,bestnest,Lb,Ub);   
%      [fnew,best,nest,fitness,bestL]=get_best_nest(nest,new_nest,fitness);
%      gbestL=bestL;
    % Update the counter
      %N_iter=N_iter+n; 
    % Discovery and randomization
      new_nest=empty_nests(nest,Lb,Ub,pa) ;
    
    % Evaluate this set of solutions
      [fnew,best,nest,fitness,Lbest]=get_best_nest(im,nest,new_nest,fitness,Lfitness);
    % Update the counter again
    %  N_iter=N_iter+n;
    % Find the best objective so far  
    if fnew<fmin,
        fmin=fnew;
        bestnest=best;
        gbestL=Lbest;
    end
    N_iter;
    fmin;
    
    fitmin(N_iter)=fmin;
    N_iter=N_iter+1;
end %% End of iterations
toc/60;


%% Post-optimization processing
%% Display all the nests
% disp(strcat('Total number of iterations=',num2str(N_iter)));
% fmin
% bestnest
% figure(12)
% plot(fitmin)
%bestnest(1,1)=400; bestnest(1,2)=40;bestnest(1,3)=2;
k=round(bestnest(1,1));
m=round(bestnest(1,2));
seRadius=bestnest(1,3);

[l,~, Sp, ~] = mslic(im, k, m, seRadius);

figure(1)
%subplot(1,2,1),
imshow(im),title('Leaf image')
%subplot(1,2,2),
show(drawregionboundaries(l, im, [255 0 255]))

end

%xlswrite('125_0019_0125_ORG.xls',fitmin');
%% --------------- All subfunctions are list below ------------------
%% Get cuckoos by ramdom walk
function nest=get_cuckoos(nest,best,Lb,Ub)
% Levy flights
n=size(nest,1);
% Levy exponent and coefficient
% For details, see equation (2.21), Page 16 (chapter 2) of the book
% X. S. Yang, Nature-Inspired Metaheuristic Algorithms, 2nd Edition, Luniver Press, (2010).
beta=3/2;
sigma=(gamma(1+beta)*sin(pi*beta/2)/(gamma((1+beta)/2)*beta*2^((beta-1)/2)))^(1/beta);

for j=1:n,
    s=nest(j,:);
    % This is a simple way of implementing Levy flights
    % For standard random walks, use step=1;
    %% Levy flights by Mantegna's algorithm
    u=randn(size(s))*sigma;
    v=randn(size(s));
    step=u./abs(v).^(1/beta);
  
    % In the next equation, the difference factor (s-best) means that 
    % when the solution is the best solution, it remains unchanged.     
    stepsize=0.01*step.*(s-best);
    % Here the factor 0.01 comes from the fact that L/100 should the typical
    % step size of walks/flights where L is the typical lenghtscale; 
    % otherwise, Levy flights may become too aggresive/efficient, 
    % which makes new solutions (even) jump out side of the design domain 
    % (and thus wasting evaluations).
    % Now the actual random walks or flights
    %s=s+stepsize.*randn(size(s));
    s(1)=s(1)+stepsize(1)*rand()*10;
    s(2)=s(2)+stepsize(2)*rand()*2;
    s(3)=s(3)+stepsize(3)*rand()*02;
   % Apply simple bounds/limits
   %nest(j,:)=s;
   nest(j,:)=simplebounds(s,Lb,Ub);
end
end

%% Find the current best nest
function [fmin,best,nest,fitness,Lbest]=get_best_nest(im,nest,newnest,fitness,Lfitness)
% Evaluating all new solutions
%global bestL

for j=1:size(nest,1),
    %fnew=fobj(newnest(j,:));
    %[fnew,l] = objective_cir_slic(newnest(j,:));
    [fnew,l] = objective_slic(im,newnest(j,:));
    if fnew<=fitness(j),
       fitness(j)=fnew;
       nest(j,:)=newnest(j,:);
       Lfitness(j,:,:)=l;
    end
end
% Find the current best
[fmin,K]=max(fitness) ;
best=nest(K,:);
Lbest(:,:)=Lfitness(K,:,:);
end

%% Replace some nests by constructing new solutions/nests
function new_nest=empty_nests(nest,Lb,Ub,pa)
% A fraction of worse nests are discovered with a probability pa
n=size(nest,1);
% Discovered or not -- a status vector
K=rand(size(nest))>pa;

% In the real world, if a cuckoo's egg is very similar to a host's eggs, then 
% this cuckoo's egg is less likely to be discovered, thus the fitness should 
% be related to the difference in solutions.  Therefore, it is a good idea 
% to do a random walk in a biased way with some random step sizes.  
%% New solution by biased/selective random walks
stepsize=rand*(nest(randperm(n),:)-nest(randperm(n),:));
new_nest=nest+stepsize.*K;
for j=1:size(new_nest,1)
    s=new_nest(j,:);
  new_nest(j,:)=simplebounds(s,Lb,Ub);  
end
end

% Application of simple constraints
function s=simplebounds(s,Lb,Ub)
  % Apply the lower bound
  ns_tmp=s;
  for k=1:3
      if(ns_tmp(k)>Ub(k))
          ns_tmp(k)=Lb(k)+(Ub(k)-Lb(k))*rand();
      end    
  
      if(ns_tmp(k)<Lb(k))
          ns_tmp(k)=Lb(k)+(Ub(k)-Lb(k))*rand();
      end 
  end
  s=ns_tmp;
%   I=ns_tmp(i)<Lb(i);
%   ns_tmp(I)=Lb(I);
%   
%   % Apply the upper bounds 
%   J=ns_tmp>Ub;
%   ns_tmp(J)=Ub(J);
%   % Update this new move 
%   s=ns_tmp;
end

%% You can replace the following by your own functions
% A d-dimensional objective function
function z=fobj(u)
%% d-dimensional sphere function sum_j=1^d (u_j-1)^2. 
%  with a minimum at (1,1, ...., 1); 
z=sum((u-1).^2);
end

function [fval,l] = objective_slic(im,pval)
%global im;
k=round(pval(1,1));
m=round(pval(1,2));
seRadius=pval(1,3);
[l,~, Sp, ~] = mslic(im, k, m, seRadius);
%show(drawregionboundaries(l, im, [255 255 255]))

% aL=Sp.stdL;
% aa=Sp.stda;
% ab=Sp.stdb;


[v,u]=size(Sp);

for i=1:u
aLab(i)=Sp(i).stdL+Sp(i).stda+Sp(i).stdb;
end

fval=sum(aLab);
end



% SLIC Simple Linear Iterative Clustering SuperPixels
%
% Implementation of Achanta, Shaji, Smith, Lucchi, Fua and Susstrunk's
% SLIC Superpixels
%
% Usage:   [l, Am, Sp, d] = slic(im, k, m, seRadius, colopt, mw)
%

% Copyright (c) 2013 Peter Kovesi
% www.peterkovesi.com/matlabfns/
% 
% Permission is hereby granted, free of charge, to any person obtaining a copy
% of this software and associated documentation files (the "Software"), to deal
% in the Software without restriction, subject to the following conditions:
% 
% The above copyright notice and this permission notice shall be included in 
% all copies or substantial portions of the Software.
%
% The Software is provided "as is", without warranty of any kind.

% Feb  2013
% July 2013 Super pixel attributes returned as a structure array




function [l, Am, Sp, d] = mslic(im, k, m, seRadius, colopt, mw, nItr, eim, We)
    
    if ~exist('colopt','var') || isempty(colopt), colopt = 'mean'; end
    if ~exist('mw','var')     || isempty(mw),         mw = 0;      end
    if ~exist('nItr','var')   || isempty(nItr),     nItr = 10;     end
    
    if exist('eim', 'var'), USEDIST = 0; else, USEDIST = 1; end
        
    MEANCENTRE = 1;
    MEDIANCENTRE = 2;
    
    if strcmp(colopt, 'mean')
        centre = MEANCENTRE;
    elseif strcmp(colopt, 'median')
        centre = MEDIANCENTRE;        
    else
        error('Invalid colour centre computation option');
    end
    
    [rows, cols, chan] = size(im);
    if chan ~= 3
        error('Image must be colour');
    end
    
    % Convert image to L*a*b* colourspace.  This gives us a colourspace that is
    % nominally perceptually uniform. This allows us to use the euclidean
    % distance between colour coordinates to measure differences between
    % colours.  Note the image becomes double after conversion.  We may want to
    % go to signed shorts to save memory.
    im = rgb2lab(im); 
    %im=double(im);
    % Apply median filtering to colour components if mw has been supplied
    % and/or non-zero
    if mw
        if length(mw) == 1
            mw(2) = mw(1);  % Use same filtering for L and chrominance
        end
        for n = 1:3
            im(:,:,n) = medfilt2(im(:,:,n), [mw(1) mw(1)]);
        end
    end
    
    % Nominal spacing between grid elements assuming hexagonal grid
    S = sqrt(rows*cols / (k * sqrt(3)/2));
    
    % Get nodes per row allowing a half column margin at one end that alternates
    % from row to row
    nodeCols = round(cols/S - 0.5);
    % Given an integer number of nodes per row recompute S
    S = cols/(nodeCols + 0.5); 

    % Get number of rows of nodes allowing 0.5 row margin top and bottom
    nodeRows = round(rows/(sqrt(3)/2*S));
    vSpacing = rows/nodeRows;

    % Recompute k
    k = nodeRows * nodeCols;
    
    % Allocate memory and initialise clusters, labels and distances.
    C = zeros(6,k);          % Cluster centre data  1:3 is mean Lab value,
                             % 4:5 is row, col of centre, 6 is No of pixels
    l = -ones(rows, cols);   % Pixel labels.
    d = inf(rows, cols);     % Pixel distances from cluster centres.
    
    % Initialise clusters on a hexagonal grid
    kk = 1;
    r = vSpacing/2;
    
    for ri = 1:nodeRows
        % Following code alternates the starting column for each row of grid
        % points to obtain a hexagonal pattern. Note S and vSpacing are kept
        % as doubles to prevent errors accumulating across the grid.
        if mod(ri,2), c = S/2; else, c = S;  end
        
        for ci = 1:nodeCols
            cc = round(c); rr = round(r);
            C(1:5, kk) = [squeeze(im(rr,cc,:)); cc; rr];
            c = c+S;
            kk = kk+1;
        end
        
        r = r+vSpacing;
    end
    
    % Now perform the clustering.  10 iterations is suggested but I suspect n
    % could be as small as 2 or even 1
    S = round(S);  % We need S to be an integer from now on
    
    for n = 1:nItr
       for kk = 1:k  % for each cluster

           % Get subimage around cluster
           rmin = max(C(5,kk)-S, 1);   rmax = min(C(5,kk)+S, rows); 
           cmin = max(C(4,kk)-S, 1);   cmax = min(C(4,kk)+S, cols); 
           subim = im(rmin:rmax, cmin:cmax, :);  
           assert(numel(subim) > 0)
           
           % Compute distances D between C(:,kk) and subimage
           if USEDIST
               D = dist(C(:, kk), subim, rmin, cmin, S, m);
           else
               D = dist2(C(:, kk), subim, rmin, cmin, S, m, eim, We);
           end

           % If any pixel distance from the cluster centre is less than its
           % previous value update its distance and label
           subd =  d(rmin:rmax, cmin:cmax);
           subl =  l(rmin:rmax, cmin:cmax);
           updateMask = D < subd;
           subd(updateMask) = D(updateMask);
           subl(updateMask) = kk;
           
           d(rmin:rmax, cmin:cmax) = subd;
           l(rmin:rmax, cmin:cmax) = subl;           
       end
       
       % Update cluster centres with mean values
       C(:) = 0;
       for r = 1:rows
           for c = 1:cols
              tmp = [im(r,c,1); im(r,c,2); im(r,c,3); c; r; 1];
              C(:, l(r,c)) = C(:, l(r,c)) + tmp;
           end
       end
       
       % Divide by number of pixels in each superpixel to get mean values
       for kk = 1:k 
           C(1:5,kk) = round(C(1:5,kk)/C(6,kk)); 
       end
       
       % Note the residual error, E, is not calculated because we are using a
       % fixed number of iterations 
    end
    
    % Cleanup small orphaned regions and 'spurs' on each region using
    % morphological opening on each labeled region.  The cleaned up regions are
    % assigned to the nearest cluster. The regions are renumbered and the
    % adjacency matrix regenerated.  This is needed because the cleanup is
    % likely to change the number of labeled regions.
    if seRadius
        [l, Am] = mcleanupregions(l, seRadius);
    else
        l = makeregionsdistinct(l);
        [l, minLabel, maxLabel] = renumberregions(l);
        Am = regionadjacency(l);    
    end

    % Recompute the final superpixel attributes and write information into
    % the Sp struct array.
    N = length(Am);
    Sp = struct('L', cell(1,N), 'a', cell(1,N), 'b', cell(1,N), ...
                'stdL', cell(1,N), 'stda', cell(1,N), 'stdb', cell(1,N), ...
                'r', cell(1,N), 'c', cell(1,N), 'N', cell(1,N));
    [X,Y] = meshgrid(1:cols, 1:rows);
    L = im(:,:,1);    
    A = im(:,:,2);    
    B = im(:,:,3);    
    for n = 1:N
        mask = l==n;
        nm = sum(mask(:));
        if centre == MEANCENTRE     
            Sp(n).L = sum(L(mask))/nm;
            Sp(n).a = sum(A(mask))/nm;
            Sp(n).b = sum(B(mask))/nm;
            
        elseif centre == MEDIANCENTRE
            Sp(n).L = median(L(mask));
            Sp(n).a = median(A(mask));
            Sp(n).b = median(B(mask));
        end
        
        Sp(n).r = sum(Y(mask))/nm;
        Sp(n).c = sum(X(mask))/nm;
        
        % Compute standard deviations of the colour components of each super
        % pixel. This can be used by code seeking to merge superpixels into
        % image segments.  Note these are calculated relative to the mean colour
        % component irrespective of the centre being calculated from the mean or
        % median colour component values.
        Sp(n).stdL = std(L(mask));
        Sp(n).stda = std(A(mask));
        Sp(n).stdb = std(B(mask));

        Sp(n).N = nm;  % Record number of pixels in superpixel too.
    end
end    
%-- dist -------------------------------------------
%
% Usage:  D = dist(C, im, r1, c1, S, m)
% 
% Arguments:   C - Cluster being considered
%             im - sub-image surrounding cluster centre
%         r1, c1 - row and column of top left corner of sub image within the
%                  overall image.
%              S - grid spacing
%              m - weighting factor between colour and spatial differences.
%
% Returns:     D - Distance image giving distance of every pixel in the
%                  subimage from the cluster centre
%
% Distance = sqrt( dc^2 + (ds/S)^2*m^2 )
% where:
% dc = sqrt(dl^2 + da^2 + db^2)  % Colour distance
% ds = sqrt(dx^2 + dy^2)         % Spatial distance
%
% m is a weighting factor representing the nominal maximum colour distance
% expected so that one can rank colour similarity relative to distance
% similarity.  try m in the range [1-40] for L*a*b* space
%
% ?? Might be worth trying the Geometric Mean instead ??
%  Distance = sqrt(dc * ds)
% but having a factor 'm' to play with is probably handy

% This code could be more efficient

function D = dist(C, im, r1, c1, S, m)

    % Squared spatial distance
    %    ds is a fixed 'image' we should be able to exploit this
    %    and use a fixed meshgrid for much of the time somehow...
    [rows, cols, chan] = size(im);
    [x,y] = meshgrid(c1:(c1+cols-1), r1:(r1+rows-1));
    x = x-C(4);  % x and y dist from cluster centre
    y = y-C(5);
    ds2 = x.^2 + y.^2;
    
    % Squared colour difference
    for n = 1:3
        im(:,:,n) = (im(:,:,n)-C(n)).^2;
    end
    dc2 = sum(im,3);
    
    D = sqrt(dc2 + ds2/S^2*m^2);
    
    
    
%--- dist2 ------------------------------------------
%
% Usage:  D = dist2(C, im, r1, c1, S, m, eim)
% 
% Arguments:   C - Cluster being considered
%             im - sub-image surrounding cluster centre
%         r1, c1 - row and column of top left corner of sub image within the
%                  overall image.
%              S - grid spacing
%              m - weighting factor between colour and spatial differences.
%            eim - Edge strength sub-image corresponding to im
%
% Returns:     D - Distance image giving distance of every pixel in the
%                  subimage from the cluster centre
%
% Distance = sqrt( dc^2 + (ds/S)^2*m^2 )
% where:
% dc = sqrt(dl^2 + da^2 + db^2)  % Colour distance
% ds = sqrt(dx^2 + dy^2)         % Spatial distance
%
% m is a weighting factor representing the nominal maximum colour distance
% expected so that one can rank colour similarity relative to distance
% similarity.  try m in the range [1-40] for L*a*b* space
%
end

function D = dist2(C, im, r1, c1, S, m, eim, We)

    % Squared spatial distance
    %    ds is a fixed 'image' we should be able to exploit this
    %    and use a fixed meshgrid for much of the time somehow...
    [rows, cols, chan] = size(im);
    [x,y] = meshgrid(c1:(c1+cols-1), r1:(r1+rows-1));
    x = x-C(4);
    y = y-C(5);
    ds2 = x.^2 + y.^2;
    
    % Squared colour difference
    for n = 1:3
        im(:,:,n) = (im(:,:,n)-C(n)).^2;
    end
    dc2 = sum(im,3);
    
    % Combine colour and spatial distance measure
    D = sqrt(dc2 + ds2/S^2*m^2);
    
    % for every pixel in the subimage call improfile to the cluster centre
    % and use the largest value as the 'edge distance'
    rCentre = C(5)-r1;   % Cluster centre coords relative to this sub-image
    cCentre = C(4)-c1;
    de = zeros(rows,cols);
    for r = 1:rows
        for c = 1:cols
            v = improfile(eim,[c cCentre], [r rCentre]);
            de(r,c) = max(v);
        end
    end

    % Combine edge distance with weight, We with total Distance.
    D = D + We * de;
    

end


% MCLEANUPREGIONS  Morphological clean up of small segments in an image of segmented regions
%
% Usage: [seg, Am] = mcleanupregions(seg, seRadius)
%
% Arguments: seg - A region segmented image, such as might be produced by a
%                  graph cut algorithm.  All pixels in each region are labeled
%                  by an integer.
%       seRadius - Structuring element radius.  This can be set to 0 in which
%                  case  the function will simply ensure all labeled regions
%                  are distinct and relabel them if necessary. 
%
% Returns:   seg - The updated segment image.
%             Am - Adjacency matrix of segments.  Am(i, j) indicates whether
%                  segments labeled i and j are connected/adjacent
%
% Typical application:
% If a graph cut or superpixel algorithm fails to converge stray segments
% can be left in the result.  This function tries to clean things up by:
% 1) Checking there is only one region for each segment label. If there is
%    more than one region they are given unique labels.
% 2) Eliminating regions below the structuring element size
%
% Note that regions labeled 0 are treated as a 'privileged' background region
% and is not processed/affected by the function.
%
% See also: REGIONADJACENCY, RENUMBERREGIONS, CLEANUPREGIONS, MAKEREGIONSDISTINCT

% Copyright (c) 2013 Peter Kovesi
% www.peterkovesi.com/matlabfns/
%
% Permission is hereby granted, free of charge, to any person obtaining a copy
% of this software and associated documentation files (the "Software"), to deal
% in the Software without restriction, subject to the following conditions:
% 
% The above copyright notice and this permission notice shall be included in 
% all copies or substantial portions of the Software.
%
% The Software is provided "as is", without warranty of any kind.
%
% March   2013 
% June    2013  Improved morphological cleanup process using distance map

function [seg, Am, mask] = mcleanupregions(seg, seRadius)
option = 2;
    % 1) Ensure every segment is distinct 
    [seg, maxlabel] = makeregionsdistinct(seg);
    
    % 2) Perform a morphological opening on each segment, subtract the opening
    % from the orignal segment to obtain regions to be reassigned to
    % neighbouring segments.
    if seRadius
        se = circularstruct(seRadius);   % Accurate and not noticeably slower
                                         % if radius is small
%       se = strel('disk', seRadius, 4);  % Use approximated disk for speed
        mask = zeros(size(seg));

        if option == 1        
            for l = 1:maxlabel
                b = seg == l;
                mask = mask | (b - imopen(b,se));
            end
            
        else   % Rather than perform a morphological opening on every
               % individual region in sequence the following finds separate
               % lists of unconnected regions and performs openings on these.
               % Typically an image can be covered with only 5 or 6 lists of
               % unconnected regions.  Seems to be about 2X speed of option
               % 1. (I was hoping for more...)
            list = finddisconnected(seg);
            
            for n = 1:length(list)
                b = zeros(size(seg));
                for m = 1:length(list{n})
                    b = b | seg == list{n}(m);
                end

                mask = mask | (b - imopen(b,se));
            end
        end
        
        % Compute distance map on inverse of mask
        [~, idx] = bwdist(~mask);
        
        % Assign a label to every pixel in the masked area using the label of
        % the closest pixel not in the mask as computed by bwdist
        seg(mask) = seg(idx(mask));
    end
    
    % 3) As some regions will have been relabled, possibly broken into several
    % parts, or absorbed into others and no longer exist we ensure all regions
    % are distinct again, and renumber the regions so that they sequentially
    % increase from 1.  We also need to reconstruct the adjacency matrix to
    % reflect the changed number of regions and their relabeling.

    seg = makeregionsdistinct(seg);
    [seg, minLabel, maxLabel] = renumberregions(seg);
    Am = regionadjacency(seg);    
    
end


% MAKEREGIONSDISTINCT Ensures labeled segments are distinct
%
% Usage: [seg, maxlabel] = makeregionsdistinct(seg, connectivity)
%
% Arguments: seg - A region segmented image, such as might be produced by a
%                  superpixel or graph cut algorithm.  All pixels in each
%                  region are labeled by an integer.
%   connectivity - Optional parameter indicating whether 4 or 8 connectedness
%                  should be used.  Defaults to 4.
%
% Returns:   seg - A labeled image where all segments are distinct.
%       maxlabel - Maximum segment label number.
%
% Typical application: A graphcut or superpixel algorithm may terminate in a few
% cases with multiple regions that have the same label.  This function
% identifies these regions and assigns a unique label to them.
%
% See also: SLIC, CLEANUPREGIONS, RENUMBERREGIONS

% Copyright (c) 2013 Peter Kovesi
% www.peterkovesi.com/matlabfns/
% 
% Permission is hereby granted, free of charge, to any person obtaining a copy
% of this software and associated documentation files (the "Software"), to deal
% in the Software without restriction, subject to the following conditions:
% 
% The above copyright notice and this permission notice shall be included in 
% all copies or substantial portions of the Software.
%
% The Software is provided "as is", without warranty of any kind.

% June 2013


function [seg, maxlabel] = makeregionsdistinct(seg, connectivity)
    
    if ~exist('connectivity', 'var'), connectivity = 4; end
    
    % Ensure every segment is distinct but do not touch segments 
    % with a label of 0
    labels = unique(seg(:))';
    maxlabel = max(labels);
    labels = setdiff(labels,0);  % Remove 0 from the label list
    
    for l = labels
        [bl,num] = bwlabel(seg==l, connectivity);  
        
        if num > 1  % We have more than one region with the same label
            for n = 2:num
                maxlabel = maxlabel+1;  % Generate a new label
                seg(bl==n) = maxlabel;  % and assign to this segment
            end
        end
    end

end


% CIRCULARSTRUCT
%
% Function to construct a circular structuring element
% for morphological operations.
%
% function strel = circularstruct(radius)
%
% Note radius can be a floating point value though the resulting
% circle will be a discrete approximation
%
% Peter Kovesi   March 2000

function strel = circularstruct(radius)

if radius < 1
  error('radius must be >= 1');
end

dia = ceil(2*radius);  % Diameter of structuring element

if mod(dia,2) == 0     % If diameter is a odd value
 dia = dia + 1;        % add 1 to generate a `centre pixel'
end

r = fix(dia/2);
[x,y] = meshgrid(-r:r);
rad = sqrt(x.^2 + y.^2);  
strel = rad <= radius;

end


% FINDDISCONNECTED find groupings of disconnected labeled regions
%
% Usage: list = finddisconnected(l)
%
% Argument:   l - A labeled image segmenting an image into regions, such as
%                 might be produced by a graph cut or superpixel algorithm.
%                 All pixels in each region are labeled by an integer.
%
% Returns: list - A cell array of lists of regions that are not
%                 connected. Typically there are 5 to 6 lists.
%
% Used by MCLEANUPREGIONS to reduce the number of morphological closing
% operations 
%
% See also: MCLEANUPREGIONS, REGIONADJACENCY

% Copyright (c) 2013 Peter Kovesi
% www.peterkovesi.com/matlabfns/
% 
% Permission is hereby granted, free of charge, to any person obtaining a copy
% of this software and associated documentation files (the "Software"), to deal
% in the Software without restriction, subject to the following conditions:
% 
% The above copyright notice and this permission notice shall be included in 
% all copies or substantial portions of the Software.
%
% The Software is provided "as is", without warranty of any kind.

% PK July 2013


function list = finddisconnected(l)
 
    debug = 0;
    [Am, Al] = regionadjacency(l);
    
    N = max(l(:));  % number of labels
    
    % Array for keeping track of visited labels
    visited = zeros(N,1);

    list = {};
    listNo = 0;
    for n = 1:N

        if ~visited(n)
            listNo = listNo + 1;
            list{listNo} = n;
            visited(n) = 1;
            
            % Find all regions not directly connected to n and not visited
            notConnected = setdiff(find(~Am(n,:)), find(visited));
            
            % For each unconnected region check that it is not already
            % connected to a region in the list. If not, add to list
            for m = notConnected
                if isempty(intersect(Al{m}, list{listNo}))
                    list{listNo} = [list{listNo} m];
                    visited(m) = 1;
                end
            end
         end % if not visited(n)
        
    end
    
    % Display each list of unconncted regions as an image
    if debug   
        for n = 1:length(list)
            
            mask = zeros(size(l));
            for m = 1:length(list{n})
                mask = mask | l == list{n}(m);
            end
            
            fprintf('list %d of %d length %d \n', n, length(list), length(list{n}))
            show(mask);
            keypause
        end
    end
end


% REGIONADJACENCY Computes adjacency matrix for image of labeled segmented regions
%
% Usage:  [Am, Al] = regionadjacency(L, connectivity)
%
% Arguments:  L - A region segmented image, such as might be produced by a
%                 graph cut or superpixel algorithm.  All pixels in each
%                 region are labeled by an integer.
%  connectivity - 8 or 4.  If not specified connectivity defaults to 8.
%
% Returns:   Am - An adjacency matrix indicating which labeled regions are
%                 adjacent to each other, that is, they share boundaries. Am
%                 is sparse to save memory.
%            Al - A cell array representing the adjacency list corresponding
%                 to Am.  Al{n} is an array of the region indices adjacent to
%                 region n.
%
% Regions with a label of 0 are not processed. They are considered to be
% 'background regions' that are not to be considered.  If you want to include
% these regions you should assign a new positive label to these areas using, say
% >> L(L==0) = max(L(:)) + 1;
%
% See also: CLEANUPREGIONS, RENUMBERREGIONS, SLIC

% Copyright (c) 2013 Peter Kovesi
% www.peterkovesi.com/matlabfns/
%
% Permission is hereby granted, free of charge, to any person obtaining a copy
% of this software and associated documentation files (the "Software"), to deal
% in the Software without restriction, subject to the following conditions:
% 
% The above copyright notice and this permission notice shall be included in 
% all copies or substantial portions of the Software.
%
% February 2013  Original version
% July     2013  Speed improvement in sparse matrix formation (4x)

function  [Am, varargout] = regionadjacency(L, connectivity)

    if ~exist('connectivity', 'var'), connectivity = 8; end
    [rows,cols] = size(L);
    
    % Identify the unique labels in the image, excluding 0 as a label.
    labels = setdiff(unique(L(:))',0);

    if isempty(labels)
        warning('There are no objects in the image')
        Am = [];
        Al = {};
        return
    end

    N = max(labels);    % Required size of adjacency matrix
    
    % Strategy:  Step through the labeled image.  For 8-connectedness inspect 
    % pixels as follows and set the appropriate entries in the adjacency
    % matrix. 
    %      x - o
    %    / | \
    %  o   o   o
    %
    % For 4-connectedness we only inspect the following pixels
    %      x - o
    %      | 
    %      o  
    %
    % Becuase the adjacency search looks 'forwards' a final OR operation is
    % performed on the adjacency matrix and its transpose to ensure
    % connectivity both ways.

    % Allocate vectors for forming row, col, value triplets used to construct
    % sparse matrix.  Forming these vectors first is faster than filling
    % entries directly into the sparse matrix
    i = zeros(rows*cols,1);  % row value
    j = zeros(rows*cols,1);  % col value
    s = zeros(rows*cols,1);  % value
    
    if connectivity == 8
        n = 1;
        for r = 1:rows-1

            % Handle pixels in 1st column
            i(n) = L(r,1); j(n) = L(r  ,2); s(n) = 1; n=n+1;
            i(n) = L(r,1); j(n) = L(r+1,1); s(n) = 1; n=n+1;
            i(n) = L(r,1); j(n) = L(r+1,2); s(n) = 1; n=n+1;
            
            % ... now the rest of the column
            for c = 2:cols-1
               i(n) = L(r,c); j(n) = L(r  ,c+1); s(n) = 1; n=n+1;
               i(n) = L(r,c); j(n) = L(r+1,c-1); s(n) = 1; n=n+1;
               i(n) = L(r,c); j(n) = L(r+1,c  ); s(n) = 1; n=n+1;
               i(n) = L(r,c); j(n) = L(r+1,c+1); s(n) = 1; n=n+1;
            end
        end
        
    elseif connectivity == 4
        n = 1;
        for r = 1:rows-1
            for c = 1:cols-1
                i(n) = L(r,c); j(n) = L(r  ,c+1); s(n) = 1; n=n+1;
                i(n) = L(r,c); j(n) = L(r+1,c  ); s(n) = 1; n=n+1;
            end
        end
    
    else
        error('Connectivity must be 4 or 8');
    end
    
    % Form the logical sparse adjacency matrix
    Am = logical(sparse(i, j, s, N, N)); 
    
    % Zero out the diagonal 
    for r = 1:N
        Am(r,r) = 0;
    end
    
    % Ensure connectivity both ways for all regions.
    Am = Am | Am';
    
    % If an adjacency list is requested...
    if nargout == 2
        Al = cell(N,1);
        for r = 1:N
            Al{r} = find(Am(r,:));
        end
        varargout{1} = Al;
    end
    
end


% RENUMBERREGIONS
%
% Usage: [nL, minLabel, maxLabel] = renumberregions(L)
%
% Argument:   L - A labeled image segmenting an image into regions, such as
%                 might be produced by a graph cut or superpixel algorithm.
%                 All pixels in each region are labeled by an integer.
%
% Returns:   nL - A relabeled version of L so that label numbers form a
%                 sequence 1:maxRegions  or 0:maxRegions-1 depending on
%                 whether L has a region labeled with 0s or not.
%      minLabel - Minimum label in the renumbered image.  This will be 0 or 1.
%      maxLabel - Maximum label in the renumbered image.
%
% Application: Segmentation algorithms can produce a labeled image with a non
% contiguous numbering of regions 1 4 6 etc. This function renumbers them into a
% contiguous sequence.  If the input image has a region labeled with 0s this
% region is treated as a privileged 'background region' and retains its 0
% labeling. The resulting image will have labels ranging over 0:maxRegions-1.
% Otherwise the image will be relabeled over the sequence 1:maxRegions
%
% See also: CLEANUPREGIONS, REGIONADJACENCY

% Copyright (c) 2010 Peter Kovesi
% www.peterkovesi.com/matlabfns/
%
% Permission is hereby granted, free of charge, to any person obtaining a copy
% of this software and associated documentation files (the "Software"), to deal
% in the Software without restriction, subject to the following conditions:
% 
% The above copyright notice and this permission notice shall be included in 
% all copies or substantial portions of the Software.
%
% October  2010
% February 2013 Return label numbering range

function [nL, minLabel, maxLabel] = renumberregions(L)

    nL = L;
    labels = unique(L(:))';  % Sorted list of unique labels
    N = length(labels);
    
    % If there is a label of 0 we ensure that we do not renumber that region
    % by removing it from the list of labels to be renumbered.
    if labels(1) == 0
        labels = labels(2:end);
        minLabel = 0;
        maxLabel = N-1;
    else
        minLabel = 1;
        maxLabel = N;
    end
    
    % Now do the relabelling
    count = 1;
    for n = labels
        nL(L==n) = count;
        count = count+1;
    end
end

% DRAWREGIONBOUNDARIES Draw boundaries of labeled regions in an image
%
% Usage: maskim = drawregionboundaries(l, im, col)
%
% Arguments:
%            l - Labeled image of regions.
%           im - Optional image to overlay the region boundaries on.
%          col - Optional colour specification. Defaults to black.  Note that
%                the colour components are specified as values 0-255.
%                For example red is [255 0 0] and white is [255 255 255].
%
% Returns: 
%       maskim - If no image has been supplied maskim is a binary mask
%                image indicating where the region boundaries are.
%                If an image has been supplied maskim is the image with the
%                region boundaries overlaid 
%
% See also: MASKIMAGE

% Copyright (c) 2013 Peter Kovesi
% www.peterkovesi.com/matlabfns/
% 
% Permission is hereby granted, free of charge, to any person obtaining a copy
% of this software and associated documentation files (the "Software"), to deal
% in the Software without restriction, subject to the following conditions:
% 
% The above copyright notice and this permission notice shall be included in 
% all copies or substantial portions of the Software.
%
% The Software is provided "as is", without warranty of any kind.

% Feb 2013

function maskim = drawregionboundaries(l, im, col)
    
    % Form the mask by applying a sobel edge detector to the labeled image,
    % thresholding and then thinning the result.
%    h = [1  0 -1
%         2  0 -2
%         1  0 -1];
    h = [-1 1];  % A simple small filter is better in this application.
                 % Small regions 1 pixel wide get missed using a Sobel
                 % operator 
    gx = filter2(h ,l);
    gy = filter2(h',l);
    maskim = (gx.^2 + gy.^2) > 0;
    maskim = bwmorph(maskim, 'thin', Inf);
    
    % Zero out any mask values that may have been set around the edge of the
    % image.
    maskim(1,:) = 0; maskim(end,:) = 0;
    maskim(:,1) = 0; maskim(:,end) = 0;
    
    % If an image has been supplied apply the mask to the image and return it 
    if exist('im', 'var') 
        if ~exist('col', 'var'), col = 0; end
        maskim = maskimage(im, maskim, col);
    end
end



% SHOW - Displays an image with the right size, colors, range, and with a title.
%
% Usage:   
%         h = show(im);
%         h = show(im, figNo, title, clim, colourmap)
%
% Arguments:  im    - Either a 2 or 3D array of pixel values or the name
%                     of an image file;
%             figNo - Optional figure number to display image in. If
%                     figNo is 0 the current figure or subplot is
%                     assumed.  The default is to create a new figure.
%             title - Optional string specifying figure title. 
%                     Defaults to the variable name of the image.
%             clim  - Optional 2-vector specifying range of values to
%                     display. Useful for dealing with outlying values im
%                     your data, or for displaying data with a diverging
%                     colour map properly. Defaults to the full data range.
%         colourmap - Optional Nx3 matrix specifying a colour map.
%                     Defaults to gray(256).
%
% Returns:    h     - Handle to the figure.  This allows you to set
%                     additional figure attributes if desired.
%
% Apart from 'im' all arguments are optional and can be supplied in any
% order.  They are inferred from their types.
%
% Where possible the image is displayed as 'TrueSize', that is, pixels on the
% screen match pixels in the image.
%
% Unless you are doing a subplot (figNo==0) the window is sized to match
% the image, leaving no border,  hence saving desktop real estate.
%
% See also: COLORCET, SHOWSURF

% Copyright (c) 2000-2018 Peter Kovesi
% peterkovesi.com
% 
% Permission is hereby granted, free of charge, to any person obtaining a copy
% of this software and associated documentation files (the "Software"), to deal
% in the Software without restriction, subject to the following conditions:
% 
% The above copyright notice and this permission notice shall be included in 
% all copies or substantial portions of the Software.
%
% The Software is provided "as is", without warranty of any kind.

% October   2000  Original version
% March     2003  Mods to alow figure name in window bar and allow for subplots.
% April     2007  Proper recording and restoring of MATLAB warning state.
% September 2008  Octave compatible
% May       2009  Reworked argument handling logic for extra flexibility
% January   2013  More Octave compatibility and proper restoring of warning state (Carnë Draug)
% March     2013  Ensure grey colour map has 256 levels rather than the
%                 default 64
% October   2017  Let imagesc apply default scaling for images < 500 pixels in size.
% April     2018  Improved argument handling.  Allow colour map and display
%                 range limits to be specified.

function h = show(im, varargin)

    Octave = exist('OCTAVE_VERSION', 'builtin') == 5; % Are we running under Octave

    s = warning('query','all');                 % Record existing warning state.
    warn_state = onCleanup (@() warning(s));    % Restore warnings at the end
    warning('off');                             % Turn off warnings that might arise if image
                                                % has to be rescaled to fit on screen

    [figNo, Title, clim, colourmap] = checkargs(varargin(:));
    
    % Check case where im is an image filename rather than image data
    if ~isnumeric(im) && ~islogical(im) 
        if isempty(Title)
            Title = im;        % Use file name for title
        end
        im = imread(im);
    elseif isempty(Title)
        Title = inputname(1);  % Use variable name of image data for title
    end

    sze = max(size(im));       % Maximum dimension of image

    if figNo > 0               % We have a valid figure number
        figure(figNo);         % Reuse or create a figure window with this number
        subplot('position',[0 0 1 1]); % Use the whole window
    elseif figNo == -1
        figNo = figure;        % Create new figure window
        subplot('position',[0 0 1 1]); % Use the whole window
    end

    if ndims(im) == 2          % Apply colour map
        if isempty(clim)
            imagesc(im);
        else
            imagesc(im, clim);
        end
        
        colormap(colourmap);  
    else
        figure(2),
        imshow(im(:,:,1:3));   % Display as RGB (ignore any alpha channel)
        title('After superpixeling')
        %imwrite(im,'10_0019_0010_SPX_B.png');
    end

    if figNo == 0              % Assume we are trying to do a subplot 
        figNo = gcf;           % Get the current figure number
        axis('image'); axis('off');
        title(Title);          % Use a title rather than rename the figure
    else
        axis('image'); axis('off');
        set(figNo,'name', ['  ' Title])

        % If not running Octave and size of image > 500 plot image at 1:1
        % resolution. Otherwise we let imagesc use its default scaling.
        if ~Octave && sze > 500
            truesize(figNo);
        end
    end

    if nargout == 1
       h = figNo;
    end

end   
%----------------------------------------------------------------------
%
% Process optional arguments. If an arguments is:
% - a scalar assume it is the figure number
% - a string, assume it is the figure title
% - a 1x2 vector assume it is the range display limits
% - a Nx3 matrix assume it is a colour map

function [figNo, Title, clim, colourmap] = checkargs(arg)
    
    % Set defaults
    figNo = -1;  % Default indicates create new figure
    Title = '';
    clim = [];
    colourmap = gray(256);

    % Overwrite defaults with anything we recognize in the arguments
    for n = 1:length(arg);
        if isscalar(arg{n})
            figNo = arg{n};
            
        elseif ischar(arg{n})
            Title = arg{n};
            
        elseif isnumeric(arg{n}) && all(size(arg{n}) == [1,2])
            clim = arg{n};
            
        elseif isnumeric(arg{n}) && size(arg{n},1) > 1 && size(arg{n},2) == 3
            colourmap = arg{n};            
            
        else
            error('Unable to process arguments');
        end
    end
end



% MASKIMAGE Apply mask to image
%
% Usage: maskedim = maskimage(im, mask, col)
%
% Arguments:    im  - Image to be masked
%             mask  - Binary masking image
%              col  - Value/colour to be applied to regions where mask == 1
%                     If im is a colour image col can be a 3-vector
%                     specifying the colour values to be applied.
%
% Returns: maskedim - The masked image
%
% See also; DRAWREGIONBOUNDARIES

% Peter Kovesi
% www.peterkovesi.com/matlabfns/
%
% Feb 2013

function maskedim = maskimage(im, mask, col)
    
    [rows,cols, chan] = size(im);
    
    % Set default colour to 0 (black)
    if ~exist('col', 'var'), col = 0; end
    
    % Ensure col has same length as image depth.
    if length(col) == 1
        col = repmat(col, [chan 1]);
    else
        assert(length(col) == chan);
    end
    
    % Perform masking
    maskedim = im;
    for n = 1:chan
        tmp = maskedim(:,:,n);
        tmp(mask) = col(n);
        maskedim(:,:,n) = tmp;
    end
end