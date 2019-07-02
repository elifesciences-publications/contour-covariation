function [DATA_ARRAY,CORR_ARRAY] = NoiseRibbonCorrelationTest(REPEATS)
    
    % Test the range of output correlations produced by the noise-mixing
    % process.
    
    %% Initialize
    
    if nargin < 1
        REPEATS = 1;
    end
    
    ALPHA = 0:0.1:1;
    ALPHA = ALPHA .^ 0.7
    %ALPHA = [0 0.15 0.3 0.35 0.4 0.45 0.5 0.55 0.6 0.65 0.7 0.85 1];
    %ALPHA = [0 0.25 0.5 0.75 1];
    
    %% Set up pairs
    [p1,p2] = ind2sub([length(ALPHA) length(ALPHA)],find(tril(ones(length(ALPHA)),-1)));
    PAIRS = [p1 p2];
    fprintf('Number of pairs: %i',length(PAIRS));

    %% Run experiment

    % Data array just counts number of times each alpha value is
    % chosen.
    DATA_ARRAY = zeros(length(ALPHA),REPEATS);
    % Corr array keeps track of the corr of every chosen image for
    % later frequency histogram binning/analysis. Not ordered. Second
    % column shows the correlation NOT chosen in each pair.
    CORR_ARRAY = zeros(size(PAIRS,1),size(PAIRS,2),REPEATS);

    ALPHA_R = zeros(length(ALPHA),length(ALPHA)-1,REPEATS);

    QUIT = 0;

    % For each repeat
    for rep = 1:REPEATS
        
        ALPHA_GO = zeros(length(ALPHA),1);

        if QUIT == 1
            break;
        end

        % Randomize trial order
        random_index = randperm( size(PAIRS,1) );

        count = 0;

        for i = random_index

            % Increment number of times these stimuli have appeared
            ALPHA_GO(PAIRS(i,1)) = ALPHA_GO(PAIRS(i,1)) + 1;
            ALPHA_GO(PAIRS(i,2)) = ALPHA_GO(PAIRS(i,2)) + 1;

            count = count + 1;
            display(count)

            if QUIT == 1
                break;
            end

            % Choose left and right image terrains
            [lTer,rTer] = GetTerrains;

            % Get images
            l_im = NoiseRibbon(lTer,ALPHA(PAIRS(i,1)));
            r_im = NoiseRibbon(rTer,ALPHA(PAIRS(i,2)));


            % Calculate linear correlations
            R = CorrelateGradients(l_im,lTer,r_im,rTer);

            % Record correlations
            ALPHA_R(PAIRS(i,1),ALPHA_GO(PAIRS(i,1)),rep) = R(1);
            ALPHA_R(PAIRS(i,2),ALPHA_GO(PAIRS(i,2)),rep) = R(2);

            if R(1)>R(2)

                DATA_ARRAY(PAIRS(i,1),rep) = DATA_ARRAY(PAIRS(i,1),rep) + 1;
                CORR_ARRAY(i,:,rep) = [R(1) R(2)];                            
                
            elseif R(2)>=R(1)
                
                DATA_ARRAY(PAIRS(i,2),rep) = DATA_ARRAY(PAIRS(i,2),rep) + 1;
                CORR_ARRAY(i,:,rep) = [R(2) R(1)];                            
                
            end

        end
    end        
    
    %% Show correlation data
    
    figure(1);
    plot(ALPHA,mean(DATA_ARRAY,2));
    xlabel('Noise alpha (0.5 gamma scaled)');
    ylabel('Trials chosen');
    axis([0 1 0 length(ALPHA)-1]);

    figure(2);
    histogram(CORR_ARRAY(:),-0.15:0.05:1); % Plot all correlations
    hold on;
    histogram(CORR_ARRAY(:,1,:),-0.15:0.05:1); % Plot chosen correlations (greater)
    hold off;
    xlabel('Contour falloff correlation (rho)');
    ylabel('Total trials chosen');
    
end

function [lTer,rTer] = GetTerrains
    
    ter = {'bo','op','ip','vo1'};
    lTer = [ter{randperm(numel(ter),1)} '_' ter{randperm(numel(ter),1)}];
    rTer = [ter{randperm(numel(ter),1)} '_' ter{randperm(numel(ter),1)}];
    while strcmp(lTer,rTer)
        rTer = [ter{randperm(numel(ter),1)} '_' ter{randperm(numel(ter),1)}];
    end 

end


function R = CorrelateGradients(l_im,lTer,r_im,rTer)
    
    % Calculate correlation for each image.
    
    % Get vectors for terrain.
    
    images = {l_im, r_im};
    terrains = {lTer, rTer};
    R = [0 0];
    
    for i = [1 2]
        im = images{i};
        ter = terrains{i};
        DATA = load(sprintf('images/terrains/masks/vectors/%s.mat',ter));
        vectors = DATA.vectors;
        vec = vectors(:); % Linearize
        pixels = find(~isnan(vec)); % List of defined vector pixel indices
        lum = nan*ones(length(pixels),1); % Initialize lum list
        for p = 1:length(pixels)
            pixel = pixels(p);
            [r,c] = ind2sub(size(vectors),pixel);
            % For each vector pixel, try getting the shifted pixel INSIDE
            try
                new_r = round(r + 1*sin(vectors(r,c)));
                new_c = round(c + -1*cos(vectors(r,c)));
                lum(p) = im(new_r,new_c);                    
            catch
                continue % Skip this pixel
            end
        end
        
        % Convert vectors to distance from [pi/2]
        f_vec = abs(pi/2 - vec(pixels));
        f_vec(f_vec>pi) = 2*pi - f_vec(f_vec>pi);
        V = f_vec(~isnan(lum)); % List of vectors with defined lum
        L = lum(~isnan(lum)); % List of defined lums
        
        % Correlate V and L
        RHO = corr([V L]);
        R(i) = -RHO(2);
        
    end
    
    
end