function mclust_metrics(timeseries_fname,firings_fname,cluster_metrics_out_fname,clip_size)

%if (nargin<1) timeseries_fname='/media/disk1/tmp/mountainlab/tmp_long_term/1518f94e6026d0f86700ab9305c3a5b3c8e1a315-mountainsort.whiten-timeseries_out.tmp'; end;
%if (nargin<2) firings_fname='/home/magland/dev/fi_ss/analyses/labs/torben/TorbenTest3/output/ms3--tetrode6/firings.mda'; end;
%if (nargin<3) cluster_metrics_out_fname='tmp.json'; end;
%if (nargin<4) clip_size=50; end;

mfilepath=fileparts(mfilename('fullpath'));
addpath([mfilepath,'/../common']);
addpath([mfilepath,'/../common/jsonlab']);

fprintf('Reading...\n');
%X=readmda_block(timeseries_fname,[1,1],[4,32556*60*10]);
X=readmda(timeseries_fname);
FF=readmda(firings_fname);

times=FF(2,:);
labels=FF(3,:);

fprintf('Extracting clips...\n');
clips=ms_extract_clips2(X,times,clip_size);

fprintf('Computing features...\n');
features=extract_features_from_clips(clips);
%figure; ms_view_clusters(features,labels);

fprintf('Computing metrics...\n');
ids=unique(labels);
[lratio,isodist]=compute_mclust_metrics(features',labels,ids);

fprintf('Writing output...\n');
OO.clusters=cell(1,length(ids));
for j=1:length(ids)
    A=struct;
    A.label=ids(j);
    A.metrics.l_ratio=lratio(j);
    A.metrics.isolation_distance=isodist(j);
    OO.clusters{j}=A;
end;
savejson('',OO,cluster_metrics_out_fname);

function [L_Ratio, IsolationDistance]=compute_mclust_metrics(Features, ClusterID, ids)

% Computes the isolation distance and the L_ratio
% Uses the matlab functions mahal and chi2cdf

% Paul Masset, CSHL, 2017

NumClusters=length(ids);
NumFeatures=size(Features,2);

IsolationDistance=zeros(NumClusters,1);
L_Ratio=zeros(NumClusters,1);


for jj=1:NumClusters
    i=ids(jj);
    
    WithinCluster=Features(ClusterID==i,:); % Features for events within cluster
    OutsideCluster=Features(ClusterID~=i,:);  % Features for events outside of cluster
    
    nWithin=sum(ClusterID==i); % Number of spikes within cluster
    nOutside=sum(ClusterID~=i); % Number of spikes outside cluster
    
    OutsideDistances=mahal(OutsideCluster,WithinCluster);  % Compute Mahalanobis distances of spikes outside that cluster
    InsideDistances=mahal(WithinCluster,WithinCluster);  % Compute Mahalanobis distances of spikes inside that cluster
    
    % Compute the L_Ratio
    L_Outside=sum(1-chi2cdf(OutsideDistances,NumFeatures)); % Compute L 
    
    L_Ratio(jj)=L_Outside/nWithin; % Normalize
    
    % Compute the isolation distance
    if nWithin<=nOutside % Only run if current cluster is less than half total spikes
        
        [SortedDistances,idxDistances]=sort(OutsideDistances); % Sort the distances
        
        IsolationDistance(jj)=SortedDistances(nWithin); % Pick the nWithin closest spikes
        
    else
        IsolationDistance(jj)=NaN; % Cant conpute if current cluster is more than half total spikes
    end
    
end

function features=extract_features_from_clips(clips)
[M,T,L]=size(clips);
npca_per_ch=2;
npca=M*npca_per_ch;
features=zeros(npca,L);
for m=1:M
    tmp=ms_event_features(clips(m,:,:),npca_per_ch);
    features((m-1)*npca_per_ch+(1:npca_per_ch),:)=tmp;
end;
