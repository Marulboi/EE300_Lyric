clc
clear all 
close all
%%
pig_num = num2str(5);
load(['C:\Users\emir.ege-nemutlu\Desktop\resp\Pig\Pig',pig_num,'\pig',pig_num,'_correctGeo_corr.mat'])
log_file = ['C:\Users\emir.ege-nemutlu\Desktop\resp\Pig\Pig',pig_num,'\SignalsCSV\png\infolog.log'];
[pace, top_vals, bot_vals] = parse_resp_log(log_file);

triHeart = TriRep(meshes.sock.faces,meshes.sock.verts);


freq = 2048;
offset=10;


for k = 1:length(pace)
    path1 = ['C:\Users\emir.ege-nemutlu\Desktop\resp\Pig\Pig',pig_num,'\Signals\', pace{k}, '.mat'];
    if length(top_vals{1}) < 2 || length(bot_vals{1}) < 2
        disp(['Skipped :', pace{k}]);
        continue
    end
    disp(['Computing :', pace{k}]);

    [folder, name, ~] = fileparts(path1);
    folder2 = strrep(folder, 'Signals', 'Signal_select_auto');
    filename2 = [name '_Signal_select_Beat.mat'];
    path2 = fullfile(folder2, filename2);
    if exist(path1, 'file') && exist(path2, 'file')
        load(path1);
        load(path2);
    else
        % One or both files do not exist — skip or continue
        fprintf('Skipping: %s or %s not found.\n', path1, path2);
        continue;  % or use 'continue' if inside a loop
    end
    load(path1); load(path2);

    if length(top_vals) > length(bot_vals)
        exhalationpeaks = top_vals{k};
        inhlationpeaks = bot_vals{k};
    else
        inhlationpeaks = top_vals{k};
        exhalationpeaks = bot_vals{k};
    end

    [timeInterval, ~, respcycle] = timeintervalgenerator2(Signal_select_Beat, inhlationpeaks, exhalationpeaks);
    [~, meanhelper] = timeintervalgenerator(Signal_select_Beat);
    N_time = length(timeInterval);

%%



    heartpot = rec.sock.Ve_filtered;


%%
ActivationPoint = zeros(3,length(Signal_select_Beat.QRS_Onset_auto));
ActivationPointInhale = zeros(3,length(inhlationpeaks));
ActivationPointExhale = zeros(3,length(exhalationpeaks));


for i = 1:length(Signal_select_Beat.QRS_Onset_auto)
    QRS = Signal_select_Beat.QRS_Onset_auto(i):Signal_select_Beat.QRS_Onset_auto(i)+Signal_select_Beat.QRS_range;
    QRST = Signal_select_Beat.QRS_Onset_auto(i):Signal_select_Beat.QRS_Onset_auto(i)+Signal_select_Beat.Beat_range;
    
    mdVdT = getTimesElectrogram(heartpot,QRS,[]);
    

    AT(i,:) = mdVdT - QRS(i);
    AT_ms(i,:) = (mdVdT - QRS(i))/(freq/1000);
    % AT_JD(i,:) = depolMapJD(triHeart,heartpot(:,QRS));
  

 
    mdVdT(rec.sock.BadChannels) = nan;
    badbyeeye = Signal_select_Beat.badATLeads;
    mdVdT(badbyeeye) = nan;

    % Signal_select_Beat.badATLeads = badbyeeye;

    interpAt = medianSpatialFilter(triHeart.X(isnan(mdVdT),:), triHeart.X(~isnan(mdVdT),:),mdVdT(~isnan(mdVdT))',10); % 30 max 
    inter_mdVdT =mdVdT;
    inter_mdVdT(isnan(mdVdT)) = interpAt;




    %% Laura's function torn apart and changed
    % Define if single point or group of points have earliest activation 
    [ATOrdered, Order] = sort(inter_mdVdT);
    [earlySites] = find(ATOrdered==ATOrdered(1));
    if(isscalar(earlySites)) %  Single early activation
        paceingsite =  triHeart.X(Order(1),:);
    else % Multiple points with same activation, try to group them
        group{1} = Order(1);
        k1=1;n=1;
        Order = Order((1:length(earlySites)));
        OrderAll = Order;
        while(~isempty(Order))
            [neighbour, distance] = findNearNeighbours(triHeart,Order(1));
            neighbour = intersect(neighbour,OrderAll);
            if(isempty(neighbour))
                % single point with no neighbours
                k1 = k1+1;
                group{k1} =  Order(1);
                Order(1)=[];
            else
                if(isempty(intersect([neighbour; Order(1)],group{k1})))
                    % check other groups also not connected
                    m = k1-1;
                    while(m>0)
                        %disp('second while loop')
                        if(isempty(intersect([neighbour; Order(1)],group{m})))
                            m = m-1;
                        else
                            group{m} = unique([group{m} neighbour']);
                            Order(1)=[];
                            m=-1;
                        end
                    end
                    if(m==0)
                        k1 = k1+1;
                        group{k1} =  Order(1);
                        Order(1)=[];
                    end
                else
                    group{k1} = unique([group{k1} neighbour']);
                    Order(1)=[];
                end
            end
        end
        % Defines pacing site as the centre of the biggest group
        if(size(group,2)==1)
            if(size(group{1},2)>1)
                paceingsite =  mean(triHeart.X(group{1},:));
            else
                paceingsite =  triHeart.X(group{1},:);
            end
        else % multiple groups; define biggest
            [~,id] = max(cellfun(@numel,group));
            if(size(group{id},2)>1)
                paceingsite =  mean(triHeart.X(group{id},:));
            else
                paceingsite =  triHeart.X(group{id},:);
            end  
        end        
    end

    ActivationPoint(:,i) = paceingsite';
    if respcycle(i) == 1
        ActivationPointInhale(:,i) = paceingsite';
    elseif respcycle(i) == 2
        ActivationPointExhale(:,i) = paceingsite';
    end
    ActivationPointInhale(:,all(ActivationPointInhale == 0))=[];
    ActivationPointExhale(:,all(ActivationPointExhale == 0))=[];

end
disp('=========================AllByTime==================')
disp(ActivationPoint)

disp('=========================Inhale==================')
disp(ActivationPointInhale)
disp('=========================Exhale==================')
disp(ActivationPointExhale)
%save(beatdata, 'Signal_select_Beat')




%% output AT maps

numPoints = size(ActivationPoint, 2);
colors = lines(numPoints); % Generate distinguishable colors


close all;
Dessin(triHeart,inter_mdVdT,inter_mdVdT)
hold on
for i1 = 1:numPoints
    scatter3(ActivationPoint(1,i1), ActivationPoint(2,i1), ActivationPoint(3,i1), ...
             100, colors(i1,:), 'filled');  % 50 is marker size
end

pause();

end