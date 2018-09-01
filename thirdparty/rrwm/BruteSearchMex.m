% 
% BruteSearchMex 
% 
%     When the dataset is small, when you have to run only a few number of
%     search, or when the dimensions of points is large, the brute search
%     method is still faster than kd-trees data structure. Computing the
%     distances one by one take a minor time than building the kd-tree. I
%     saw many k-neighbours utilities on FEX, but all them were m-coded. I
%     think such a brute calculation is not an m-code job. So I developed
%     my own Nearest Neighbour finder. It is nothing special it just
%     computes all the distances and take the ones required form the input
%     parameters, but of course the mex implementation make it even faster
%     than vectorized m-code. 
% 
%     The following utilities are provided:
%     
%     - Nearest neighbor
%     
%     - K-Nearest neighbors
%     
%     - Radius Search
%     
%     
%     They al supports N-dimensions and work on double, it is possible to
%     choose if return the distances.
% 
% Here are some example that explain all possible sintax:
% 
%N=10000;%number of reference points
%Nq=1000;%number of query points
%dim=10;%dimension of points
%k=3;%number of neighbors
%r=.1;%Search radius
% 
% p=rand(N,dim);%reference points 
% qp=rand(Nq,dim);%query points
% 
% 
% %To choose the search method is just necessary to add a string and a
% %parameter
% 
%% Nearest neighbor
% 
% idc=BruteSearchMex(p,qp);%find the nearest neighbour for each query
% points
% 
% [idc,dist]=BruteSearchMex(p,qp);%same but returns also the distance of
% pints
% 
%
%
%% K-Nearest neighbor
% 
% kidc=BruteSearchMex(p,qp,'k',k);%find the K nearest neighbour for each
% query points
% 
% [kidc,kdist]=BruteSearchMex(p,qp,'k',k);%same but returns also the distance
% of pints
% 
% 
%% Radius Search
% 
% %NOTE: Differently from the others the radius search only supports one
% input query %point
%for i=1:Nq
% 
%     ridc=BruteSearchMex(p,qp(i,:),'r',r);%find thepoints within the
%     distance of the radius from query point
% 
%     [ridc,rdist]=BruteSearchMex(p,qp(i,:),'r',r);%same but returns also the
%     distance of pints
% 
% end
%
% See also, kdtree, nnsearch, delaunary, dsearch
%
% For info, questions, suggestions, bugs: giaccariluigi@msn.com
%
%
%Author: Giaccari Luigi.
%Last Update 10/12/2008.
%Created:    10/12/2008.
