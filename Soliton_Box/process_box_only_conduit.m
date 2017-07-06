%% Extract box peaks; get number, ampl dist
load('box_num_results.mat','ps','ampls','ws','asolis','boxpks');
numruns = numel(ws);
N_num   = zeros(1,numruns);
ampls_num = struct;

for ii = 1:numruns
    b = ['run',num2str(ii)];
    bm = [];
   % Turn boxpks into a matrix
   for jj = 1:10
       [m,n] = size(bm);
       q = length(boxpks.(b).(['t',num2str(jj)]));
       if m<q
            bm = [bm; zeros(q-m,n)];
       elseif m>q
           boxpks.(b).(['t',num2str(jj)]) = [boxpks.(b).(['t',num2str(jj)]); zeros(m-q,1) ];
       end
        bm = [bm, boxpks.(b).(['t',num2str(jj)]) ];
   end
    [m,n] = size(bm);
    N_num(ii) = m;
    ampls_num.(b) = mean(bm,2);
end

save('box_num_results.mat','N_num','ampls_num','-append');