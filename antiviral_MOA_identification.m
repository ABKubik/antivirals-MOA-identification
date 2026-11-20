%% Fitting the data from POSITIVES
% (observational cohort study Post-vaccination Viral Characteristics Study)

% This code generates the results presented in "Discovering the modes of
% action of antivirals from data and modelling: SARS-CoV-2 and Paxlovid as
% case study", S.M. Ciupe, A.B. Kubik, G. Röst (2026)

% All auxiliary functions can be found at the end of this code

clear all, close all

% In the following line, one should load the de-identified viral load data
% that can be found in https://github.com/tinphan-LANL. This repository is
% reported in "Modeling suggests SARS-CoV-2 rebound after nirmatrelvir-
% ritonavir treatment is driven by target cell preservation coupled with
% incomplete viral clearance", T. Phan, R.M. Ribeiro, G.E. Edelstein et al.
% (Jorunal of Virology, 99(3), 2025)

% T = readtable('...');

% We extract the data of Groups 1, 2 and 4 from the table
T1 = T(T.Group == 1, :); % Treated group, rebound
T2 = T(T.Group == 2, :); % Treated group, no rebound
T4 = T(T.Group == 4, :); % Untreated group, no rebound
T = [T1; T2; T4]; % New table we work with

% We collect the IDs of the patients
PIDs = unique(T.PID);
num_IDs = length(PIDs); % Number of patients

% We define the minimum and maximum days for which we have data, and define
% the time span
min_day = min(0, min(T.DaysPostSymptomOnset));
max_day = max(T.DaysPostSymptomOnset);

t_span = min_day:max_day;

% We save the curves of Treated Group depending on the day of starting
% treatment, {0, 1, 2, 3, 4, 5-6-9}, and the curves of Untreated Group
Days_vs_VL = NaN*ones(10, max_day - min_day + 1, num_IDs);

for i = 1:num_IDs

    patientData = T(T.PID == PIDs(i), :);

    if patientData.Group == 1

        measuredDays = patientData.DaysPostSymptomOnset - min_day + 1; % Adjust for indexing
        viralLoads = patientData.Log10VL; % Extract viral loads
        day_NR = patientData.DayN_RStartedAfterSymptomOnset(1) + 1; % Adjust for indexing

        Days_vs_VL(day_NR, measuredDays, i) = viralLoads; % Store viral loads

    elseif patientData.Group == 2

        measuredDays = patientData.DaysPostSymptomOnset - min_day + 1; % Adjust for indexing
        viralLoads = patientData.Log10VL; % Extract viral loads
        day_NR = patientData.DayN_RStartedAfterSymptomOnset(1) + 1; % Adjust for indexing

        Days_vs_VL(min(day_NR,6) + 3, measuredDays, i) = viralLoads; % Store viral loads

    else

        measuredDays = patientData.DaysPostSymptomOnset - min_day + 1; % Adjust for indexing
        viralLoads = patientData.Log10VL; % Extract viral loads
        % From first viral load = 1 (on log10), we disregard the next data
        inds1 = find(viralLoads==1);

        if ~isempty(inds1)

            first1ind = inds1(1);
            viralLoads(first1ind+1:end) = NaN;

        end
        
        Days_vs_VL(10, measuredDays, i) = viralLoads; % Store viral loads

    end

end


%% Fitting Least Squares for Untreated Group
% We assume knowledge of (T(-3),I(-3),V(-3))=(4e+8,0,0.31), as in "In-host
% mathematical modelling of COVID-19 in humans", E.A. Hernandez-Vargas,
% J.X. Velasco-Hernandez (Annual Reviews in Control, 50:448-456, 2020)

T0 = 4e+8; I0 = 0; V0 = 0.31;

% Time step for Runge-Kutta 4 integrator
hstep = 2^(-9);

% Chosen lower and upper bounds for (beta, delta, p, c)
lower_bound = [1e-7, 0.61, 0.1, 1];
upper_bound = [1e-6, 2.01, 1, 15]; 

PIDsG4 = unique(T.PID(T.Group == 4));
group4_ind = find(PIDs == PIDsG4(1)); % First ID of Untreated Group
max_day4 = 30; % Maximum day of simulation

len = length(-3:hstep:max_day4); 
sum_Group4UG = zeros(len, 1); % Adding up to compute mean of curves afterwards

% Table to save the fitting parameters, error and AIC of Untreated Group
dataTable = table('Size', [num_IDs-group4_ind+2, 8], ...
    'VariableTypes', ["double", "double", "double", "double", "double", ...
    "double", "double", "double"],...
    'VariableNames', ["pat_ID", "beta", "delta", "p", "c", "R0",...
    "Error", "AIC"]);

% Table with suitable initial seeds to initialize the fitting algorithm
IS_table = readtable('initialSeeds.xlsx');

% TIV system
TIV = @(t,x,params) [- params(1)*x(1)*x(3);
    params(1)*x(1)*x(3) - params(2)*x(2);
    params(3)*x(2) - params(4)*x(3)];

peak_values = zeros(1,42); % We save peak viral values...
peak_times = zeros(1,42); % ... as well as their times
clearance_times = zeros(1,42); % We save virus clearance times

% We go through the IDs of Group 4
for i = group4_ind:num_IDs
       
    figure(2)
    % 6x7 subfigures for 41 patients
    subplot(6,7,mod(i-9,42))
    hold on

    % disp(i)
    
    dataTable(i-group4_ind+1, 1) = {PIDs(i)}; % Patient ID

    % We collect and plot patient data
    patientData = T(T.PID == PIDs(i), :);
        
    day_NR = -1; % To indicate no treatment
    patientCurve = Days_vs_VL(10, :, i); % We recover the data we saved at the beginning
    
    % We collect non-NaN data
    t_span_scarce = t_span(~isnan(patientCurve));
    patientCurve_scarce = patientCurve(~isnan(patientCurve));
    n = length(patientCurve_scarce);

    t_span_Mes = t_span_scarce(patientCurve_scarce ~= 1); % Measurable data > 1
    patientCurve_Mes = patientCurve_scarce(patientCurve_scarce ~= 1);
    t_span_nonMes = t_span_scarce(patientCurve_scarce == 1); % Unmeasurable data == 1

    plot(t_span_Mes, patientCurve_Mes, 'ko', 'MarkerFaceColor', 'k',...
        'HandleVisibility','off')
    plot(t_span_nonMes, patientCurve_scarce(ismember(t_span_scarce, t_span_nonMes)),...
        'ko','HandleVisibility','off')
    plot([-5, max_day4+5], [1, 1], '--', 'Color', [0.5, 0.5, 0.5], ...
        'HandleVisibility', 'off')
    axis([-5, max_day4+5, -5, 10])

    % This function gives the error of the fitting given some parameters
    fun_norm = @(params) error_fun2([-3, t_span_scarce(end)], t_span_Mes, ...
        t_span_nonMes, patientCurve_Mes, [T0, I0, V0], hstep, params);

    % We take the initial condition from table IS_table
    IS = table2array(IS_table(i - group4_ind + 1, 2:5));

    % We compute optimal fitting:
    % bestParsUG: best parameters
    % final_errTV: final error of the fitting
    [bestParsUG, final_errUG] = fmincon(fun_norm, IS, [], [], [], [], ...
        lower_bound, upper_bound);
    
    dataTable(i-group4_ind+1, 2) = {bestParsUG(1)}; % beta
    dataTable(i-group4_ind+1, 3) = {bestParsUG(2)}; % delta
    dataTable(i-group4_ind+1, 4) = {bestParsUG(3)}; % p
    dataTable(i-group4_ind+1, 5) = {bestParsUG(4)}; % c
    dataTable(i-group4_ind+1, 6) = {bestParsUG(1)*bestParsUG(3)*T0/(bestParsUG(4)*bestParsUG(2))}; % R0
    dataTable(i-group4_ind+1, 7) = {final_errUG}; % error
    dataTable(i-group4_ind+1, 8) = {n*log(final_errUG^2/n) + 2*4}; % AIC
    % AIC = n*ln(error^2/n) + 2k, where n = #data and k = #parameters

    % We compute the solution associated to the optimal fitting and plot it
    [tMVcurve, xUG] = rungeKutta4(@(t,x) TIV(t, x, bestParsUG), ...
        [-3, max_day4], [T0, I0, V0], hstep);
    logxUG = log10(xUG(:, 3)); % log10 to compare with reported data

    plot(tMVcurve, logxUG, 'b-.', 'LineWidth', 1.5, 'DisplayName', ...
        'Individual fit')
    
    % We save the other values of interest: peak value and time, and
    % clearance time
    peak_values(i-group4_ind+1) = max(logxUG);
    peak_times(i-group4_ind+1) = tMVcurve(logxUG == max(logxUG));

    clearance_inds_log = find(logxUG-1 <= 0);
    clearance_inds_time = find(tMVcurve >= 5); % We need to check later times
    % since viral load can be below 1 also at the beginning of the infection
    clearance_1ind = find(ismember(clearance_inds_log, clearance_inds_time), 1);
    clearance_times(i-group4_ind+1) = tMVcurve(clearance_inds_log(clearance_1ind));

    sum_Group4UG = sum_Group4UG + logxUG; % We add up all fitted curves
    % in order to compute the mean curve afterwards

    s = sprintf('PID %d', PIDs(i));
    title(s)

end


%% Mean of Untreated Group curves
% We compute the mean of the fitted curves and fit it as well

dataTable(num_IDs-group4_ind+2, 1) = {97}; % Artificial ID for saving
% mean curve data in the table

% Computing the mean curve
mean_Group4UG = sum_Group4UG/(num_IDs-group4_ind+1);
n = length(mean_Group4UG);

% Fitting the mean curve
fun_meanUG = @(params) norm(log10(rungeKutta4aux(@(t,x) TIV(t,x,params),...
    [-3, max_day4], [T0, I0, V0], hstep)) - mean_Group4UG);

% We take the initial condition from table IS_table
IS = table2array(IS_table(num_IDs - group4_ind + 2, 2:5));

[fitMean, errMean] = fmincon(fun_meanUG, IS, [], [], [], [], lower_bound, ...
    upper_bound);

dataTable(num_IDs-group4_ind+2, 2) = {fitMean(1)}; % beta
dataTable(num_IDs-group4_ind+2, 3) = {fitMean(2)}; % delta
dataTable(num_IDs-group4_ind+2, 4) = {fitMean(3)}; % p
dataTable(num_IDs-group4_ind+2, 5) = {fitMean(4)}; % c
dataTable(num_IDs-group4_ind+2, 6) = {fitMean(1)*fitMean(3)*T0/(fitMean(4)*fitMean(2))}; % R0
dataTable(num_IDs-group4_ind+2, 7) = {errMean}; % error
dataTable(num_IDs-group4_ind+2, 8) = {n*log(errMean^2/n)+2*4}; % AIC

[tMVcurve, MVcurve] = rungeKutta4(@(t,x) TIV(t,x,fitMean), [-3, max_day4], ...
    [T0, I0, V0], hstep);
logMVcurve = log10(MVcurve(:, 3));

% We save the other values of interest: peak value and time, and
% clearance time
peak_values(end) = max(logMVcurve);
peak_times(end) = tMVcurve(logMVcurve == max(logMVcurve));
clearance_inds_log = find(logMVcurve-1 <= 0);
clearance_inds_time = find(tMVcurve >= 5);
clearance_1ind = find(ismember(clearance_inds_log,clearance_inds_time), 1);
clearance_times(end) = tMVcurve(clearance_inds_log(clearance_1ind));

% We plot the obtained mean curve in the individual plot of each patient
for i = group4_ind:num_IDs

    figure(2)

    subplot(6,7,mod(i-9,42))
    plot(tMVcurve, logMVcurve, '-', 'Color', [0.96,0.47,0.16], 'LineWidth', 1.5, ...
        'DisplayName', 'Fit of mean')

end

legend

% Create table 
filename = sprintf('UntreatedGroupFittings.xlsx');
writetable(dataTable, filename);


%% One mode of action fittings initialized with Untreated Group parameters

% Table to save the fitting parameters, delay, fitting error and AIC for
% both one MOA case and two MOAs case
dataTableG12 = table('Size', [112, 8], 'VariableTypes',...
    ["double", "double", "double", "double", "double", "double",...
    "double", "double"], 'VariableNames', ["pat_ID", "beta", "delta",...
    "p", "c", "tau", "Error", "AIC"]);

% We choose the IDs in Treated Group that have at least 3 data points or 2
% data points + 1 unmeasured data between (start of treatment) and (end of
% treatment + 3 days), and there is no rebound in these data
TG_IDs = [3, 5, 7, 22, 23, 24, 25, 27, 29, 34, 36, 39, 46, 47, 49, 50];

len12 = length(TG_IDs);

group2_ind = find(PIDs == 21); % We find the first ID of Group 2 in order
% to collect data from Days_vs_VL

% Colors for the plots when fitting {beta, delta, p, c}
colorsPars = {'red', 'blue', [0.23,0.67,0.2], [0.93,0.69,0.13]};

daysAfterTreat = 3; % We consider the treatment stops acting after 3 days 
% after last dose of medication

% Boundaries to fit each parameter {beta, delta, p, c}:
% - p and beta must decrease, so the upper bounds are the fitted values for
% curve MV
% - c and delta must increase, so the lower bounds are the fitted values
% for curve MV
low_bnds = [1e-10, fitMean(2), 1e-3, fitMean(4)];
up_bnds = [fitMean(1), 3, fitMean(3), 24];

j = 0; % Counter for figures

% We save the best values of delta and p for posterior reporting
best_delta1 = []; bd1 = 1;
best_p1 = []; bp1 = 1;

% We save the best values of tau for posterior reporting
best_tau1 = []; bt1 = 1;

% We go through the IDs of Treated Group
for i = 1:group4_ind-1
    
    % We check if the ID is in the above list TG_IDs
    if ismember(PIDs(i),TG_IDs)
        
        j = j+1;
    
        % disp(i)
            
        % Collect patient data and plot it
        patientData = T(T.PID == PIDs(i), :);
        
        day_NR = patientData.DayN_RStartedAfterSymptomOnset(1) + 1; % Day NR started
        patientCurve = Days_vs_VL((PIDs(i)<group2_ind)*day_NR +...
            (PIDs(i)>=group2_ind)*(min(day_NR,6)+3), :, i); % Collect data previously saved
    
        % Select non-NaN data
        t_span_scarce = t_span(~isnan(patientCurve));
        patientCurve_scarce = patientCurve(~isnan(patientCurve));

        % All data since start of symptoms until end of effect
        t_span_scarce_all = t_span_scarce(ismember(t_span_scarce,1:day_NR+3+daysAfterTreat));
        patientCurve_scarce_all = patientCurve_scarce(ismember(t_span_scarce,1:day_NR+3+daysAfterTreat));

        % Data since start of treatment until end of effect
        t_span_scarce = t_span_scarce(ismember(t_span_scarce,day_NR-1:day_NR+3+daysAfterTreat));
        patientCurve_scarce = patientCurve_scarce(ismember(t_span_scarce,day_NR-1:day_NR+3+daysAfterTreat));
        n = length(patientCurve_scarce); % Number of data available for fitting

        % We separate data between measurable (> 1) and unmeasurable (== 1) data
        t_span_Mes_all = t_span_scarce_all(patientCurve_scarce_all ~= 1);
        patientCurve_Mes_all = patientCurve_scarce_all(patientCurve_scarce_all ~= 1);
        t_span_nonMes_all = t_span_scarce_all(patientCurve_scarce_all == 1);
    
        t_span_Mes = t_span_scarce(patientCurve_scarce ~= 1);
        patientCurve_Mes = patientCurve_scarce(patientCurve_scarce ~= 1);
        t_span_nonMes = t_span_scarce(patientCurve_scarce == 1);
    
        s = sprintf('PID %d', PIDs(i));
        
        % Figure 3 plots best fittings
        figure(3)
        % 4x4 subfigures
        subplot(4,4,mod(j-1,16)+1)
        hold on
    
        plot(t_span_Mes_all, patientCurve_Mes_all, 'ko', 'MarkerFaceColor',...
            'k', 'HandleVisibility', 'off')
        plot(t_span_nonMes_all, patientCurve_scarce_all(ismember(t_span_scarce_all,t_span_nonMes_all)),...
            'ko', 'HandleVisibility', 'off')
        xregion(day_NR-1, day_NR+3, 'FaceColor', 'cyan', 'HandleVisibility', 'off')
        plot([-5, 15], [1, 1], '--', 'Color', [0.5,0.5,0.5], 'HandleVisibility', 'off')
        axis([-5, 15, 0, 10])

        title(s)
        
        % Figure 4 plots all fittings
        figure(4)

        subplot(4,4,mod(j-1,16)+1)
        hold on
    
        plot(t_span_Mes_all, patientCurve_Mes_all, 'ko', 'MarkerFaceColor',...
            'k', 'HandleVisibility', 'off')
        plot(t_span_nonMes_all, patientCurve_scarce_all(ismember(t_span_scarce_all,t_span_nonMes_all)),...
            'ko', 'HandleVisibility', 'off')
        xregion(day_NR-1, day_NR+3, 'FaceColor', 'cyan', 'HandleVisibility', 'off')
        plot([-5, 15], [1, 1], '--', 'Color', [0.5,0.5,0.5], 'HandleVisibility', 'off')
        axis([-5, 15, 0, 10])

        title(s)

        bestPars = zeros(1,4); % Best fittings quantities
        allPars = zeros(3,4); % All fittings quantities

        % We fit pi in {beta, delta, p, c} and delay tau
        for m = 1:4
            
            % We initialize the table with the basal parameters of MV curve
            dataTableG12(4*(j-1) + 1 + (m-1), 1) = {PIDs(i)};
            dataTableG12(4*(j-1) + 1 + (m-1), 2) = {fitMean(1)};
            dataTableG12(4*(j-1) + 1 + (m-1), 3) = {fitMean(2)};
            dataTableG12(4*(j-1) + 1 + (m-1), 4) = {fitMean(3)};
            dataTableG12(4*(j-1) + 1 + (m-1), 5) = {fitMean(4)};
        
            % Error of the fitting given some parameters
            err_funTVmean =@(params) treatment(fitMean, TIV, params,...
                m, [-3, day_NR+3+daysAfterTreat], hstep, patientCurve_Mes,...
                t_span_Mes, t_span_nonMes, day_NR-1, MVcurve);
            
            % Select a suitable initial seed for the fitting algorithm
            IS = table2array(IS_table(4*(j-1)+m+num_IDs-group4_ind+2, [m+1,6]));

            % Here we fit the parameter and the delay, where tau is in [0,2]
            [G12_pars, G12_err] = fmincon(err_funTVmean, IS, [], [], [],...
                [], [low_bnds(m), 0], [up_bnds(m), 2]);

            G12_AIC = n*log(G12_err^2/n)+2*length(G12_pars); % AIC

            dataTableG12(4*(j-1) + 1 + (m-1), m+1) = {G12_pars(1)}; % Parameter m
            dataTableG12(4*(j-1) + 1 + (m-1), 6) = {floor(4*G12_pars(2))/4}; % Delay, multiple of 0.25
            dataTableG12(4*(j-1) + 1 + (m-1), 7) = {G12_err}; % Error
            dataTableG12(4*(j-1) + 1 + (m-1), 8) = {G12_AIC}; % AIC

            if m==1
                % If m=1, beta is the best parameter up to now since it is
                % the only fitted parameter, so we save it in the best
                % parameter vector
    
                bestPars = [m, G12_pars, G12_AIC]; % Parameter, best fitting, AIC

            else
                % If m>1, we check whether this new fitting gives a better
                % AIC or not. If the new AIC is better, we change the best
                % parameter vector

                % We as well save the results in the all parameters matrix

                if G12_AIC < bestPars(4)
                    
                    allPars(m-1,:) = bestPars;
                    bestPars = [m, G12_pars, G12_AIC];

                else

                    allPars(m-1,:) = [m, G12_pars, G12_AIC];

                end

            end

        end

        % For posterior analysis of delay values
        best_tau1(bt1) = floor(4*bestPars(3))/4;
        bt1 = bt1+1;
        
        % For posterior analysis of best delta and p values
        if bestPars(1) == 2

            best_delta1(bd1) = bestPars(2);
            bd1 = bd1+1;

        elseif bestPars(1) == 3

            best_p1(bp1) = bestPars(2);
            bp1 = bp1+1;

        end

        % Solution with optimal parameters, initialized with curve MV
        [tTG, xTG] = optimCurves(fitMean, TIV, ...
            [bestPars(2), floor(4*bestPars(3))/4], bestPars(1), ...
            [-3, day_NR+3+daysAfterTreat], hstep, day_NR-1, MVcurve);

        logxTG = log10(xTG(:,3)); % log10 viral load
        
        % We plot the best fitting in Figures 3 and 4
        figure(3)
        subplot(4, 4, mod(j-1,16)+1)

        plot(tTG, logxTG, '-', 'Color', colorsPars{bestPars(1)}, ...
            'LineWidth', 1.5)
            
        figure(4)
        subplot(4, 4, mod(j-1,16)+1)

        plot(tTG, logxTG, '-', 'Color', colorsPars{bestPars(1)}, ...
            'LineWidth', 1.5)
       
        for m = 1:3

            % Solution with non-optimal parameters, initialized with curve MV
            [tTG, xTG] = optimCurves(fitMean, TIV,...
                [allPars(m,2), floor(4*allPars(m,3))/4], allPars(m,1),...
                [-3, day_NR+3+daysAfterTreat], hstep, day_NR-1, MVcurve);

            logxTG = log10(xTG(:,3));
            
            % If AICs between best parameter and this parameter are
            % comparable (i.e., distance <= 2), then we also consider the
            % current parameter as a best fitting
            if abs(allPars(m,4) - bestPars(4)) <= 2   

                figure(3)
                subplot(4, 4, mod(j-1,16)+1)

                plot(tTG, logxTG, '-', 'Color', colorsPars{allPars(m,1)}, ...
                    'LineWidth',1.5)

                best_tau1(bt1) = floor(4*allPars(m,3))/4;
                bt1 = bt1+1;

                if allPars(m,1) == 2

                    best_delta1(bd1) = allPars(m,2);
                    bd1 = bd1+1;

                elseif allPars(m,1) == 3

                    best_p1(bp1) = allPars(m,2);
                    bp1 = bp1+1;

                end

            end
            
            figure(4)
            subplot(4, 4, mod(j-1,16)+1)

            plot(tTG, logxTG, '-', 'Color', colorsPars{allPars(m,1)}, ...
                'LineWidth', 1.5)

        end
    
    end

end 
       
%% Two modes of action fittings initialized with Untreated Group parameters

parsVec = [1 2 4]; % Corresponding to positions of beta, delta and c in our
% regular parameters vector (beta, delta, p, c)

j = 0;

% As before, we save best delta, p and tau values for posterior analysis
best_delta2 = []; bd2 = 1;
best_p2 = []; bp2 = 1;
best_tau2 = []; bt2 = 1;

% We go through the IDs of Treated Group
for i = 1:group4_ind-1
    
    % We check if the ID is in the above list TG_IDs
    if ismember(PIDs(i),TG_IDs)
        
        j = j+1;
    
        % disp(i)
        
        % Collect patient data and plot it
        patientData = T(T.PID == PIDs(i), :);
        
        day_NR = patientData.DayN_RStartedAfterSymptomOnset(1) + 1; % Day NR started
        patientCurve = Days_vs_VL((PIDs(i)<group2_ind)*day_NR +...
            (PIDs(i)>=group2_ind)*(min(day_NR,6)+3), :, i); % Collect data previously saved
    
        % Select non-NaN data
        t_span_scarce = t_span(~isnan(patientCurve));
        patientCurve_scarce = patientCurve(~isnan(patientCurve));

        % All data since start of symptoms until end of effect
        t_span_scarce_all = t_span_scarce(ismember(t_span_scarce,1:day_NR+3+daysAfterTreat));
        patientCurve_scarce_all = patientCurve_scarce(ismember(t_span_scarce,1:day_NR+3+daysAfterTreat));

        % Data since start of treatment until end of effect
        t_span_scarce = t_span_scarce(ismember(t_span_scarce,day_NR-1:day_NR+3+daysAfterTreat));
        patientCurve_scarce = patientCurve_scarce(ismember(t_span_scarce,day_NR-1:day_NR+3+daysAfterTreat));
        n = length(patientCurve_scarce);

        % We separate data between measurable (> 1) and unmeasurable (== 1) data
        t_span_Mes_all = t_span_scarce_all(patientCurve_scarce_all ~= 1);
        patientCurve_Mes_all = patientCurve_scarce_all(patientCurve_scarce_all ~= 1);
        t_span_nonMes_all = t_span_scarce_all(patientCurve_scarce_all == 1);
    
        t_span_Mes = t_span_scarce(patientCurve_scarce ~= 1);
        patientCurve_Mes = patientCurve_scarce(patientCurve_scarce ~= 1);
        t_span_nonMes = t_span_scarce(patientCurve_scarce == 1);
    
        s = sprintf('PID %d', PIDs(i));

        % Figure 5 plots best fittings
        figure(5)
        % 4x4 subfigures
        subplot(4,4,mod(j-1,16)+1)
        hold on
    
        plot(t_span_Mes_all, patientCurve_Mes_all, 'ko', 'MarkerFaceColor',...
            'k', 'HandleVisibility', 'off')
        plot(t_span_nonMes_all, patientCurve_scarce_all(ismember(t_span_scarce_all,t_span_nonMes_all)),...
            'ko', 'HandleVisibility', 'off')
        xregion(day_NR-1, day_NR+3, 'FaceColor', 'cyan', 'HandleVisibility', 'off')
        plot([-5, 15], [1, 1], '--', 'Color', [0.5,0.5,0.5], 'HandleVisibility', 'off')
        axis([-5, 15, 0, 10])

        title(s)
        
        % Figure 6 plots all fittings
        figure(6)

        subplot(4,4,mod(j-1,16)+1)
        hold on
    
        plot(t_span_Mes_all, patientCurve_Mes_all, 'ko', 'MarkerFaceColor',...
            'k', 'HandleVisibility', 'off')
        plot(t_span_nonMes_all, patientCurve_scarce_all(ismember(t_span_scarce_all,t_span_nonMes_all)),...
            'ko', 'HandleVisibility', 'off')
        xregion(day_NR-1, day_NR+3, 'FaceColor', 'cyan', 'HandleVisibility', 'off')
        plot([-5, 15], [1, 1], '--', 'Color', [0.5,0.5,0.5], 'HandleVisibility', 'off')
        axis([-5, 15, 0, 10])

        title(s)

        bestPars = zeros(1,7); % Best fittings quantities
        allPars = zeros(2,7); % All fittings quantities

        % We fit p, pi in {beta, delta, c} and delay tau
        for m = 1:3
            
            % We initialize the table with the basal parameters of MV curve
            dataTableG12(64 + 3*(j-1) + 1 + (m-1), 1) = {PIDs(i)};
            dataTableG12(64 + 3*(j-1) + 1 + (m-1), 2) = {fitMean(1)};
            dataTableG12(64 + 3*(j-1) + 1 + (m-1), 3) = {fitMean(2)};
            dataTableG12(64 + 3*(j-1) + 1 + (m-1), 5) = {fitMean(4)};
            
            % Tuple of positions of parameters, where beta==1, delta==2,
            % p==3 and c==4
            posParams = sort([parsVec(m), 3]);

            % Error of the fitting given some parameters
            err_funTVmean =@(params) treatment(fitMean, TIV, params, ...
                posParams, [-3, day_NR+3+daysAfterTreat], hstep, ...
                patientCurve_Mes, t_span_Mes, t_span_nonMes, day_NR-1, MVcurve);
            
            % Select a suitable initial seed for the fitting algorithm
            %4*(j-1)+(m-1)+num_IDs-group4_ind+2
            %+4*(len12-1)
            IS = table2array(IS_table(3*(j-1)+m+4*len12+num_IDs-group4_ind+2, [posParams+1,6]));

            % Here we fit both parameters and the delay, where tau is in [0,2]
            [G12_pars, G12_err] = fmincon(err_funTVmean, IS, [], [], [], ...
                [], [low_bnds(posParams), 0], [up_bnds(posParams), 2]);

            G12_AIC = n*log(G12_err^2/n)+2*length(G12_pars); % AIC
            
            % We save the first parameter, following the ordering (beta,
            % delta, p, c)
            dataTableG12(64 + 3*(j-1) + 1 + (m-1), parsVec(m)+1) = {G12_pars((parsVec(m)<3) + 2*(parsVec(m)>3))};
            % We save the second parameter, following the ordering (beta,
            % delta, p, c)
            dataTableG12(64 + 3*(j-1) + 1 + (m-1), 4) = {G12_pars((parsVec(m)>3) + 2*(parsVec(m)<3))};
            dataTableG12(64 + 3*(j-1) + 1 + (m-1), 6) = {floor(4*G12_pars(3))/4}; % Delay, multiple of 0.25
            dataTableG12(64 + 3*(j-1) + 1 + (m-1), 7) = {G12_err}; % Error
            dataTableG12(64 + 3*(j-1) + 1 + (m-1), 8) = {G12_AIC}; % AIC

            if m==1
                % If m=1, (beta, p) is the best parameter tuple up to now 
                % since it is the only fitted pair, so we save it in the
                % best parameter vector
    
                bestPars = [posParams, G12_pars, G12_AIC, parsVec(m)];
                % Parameters, best fitting, AIC, parameter pi

            else
                % If m>1, we check whether this new fitting gives a better
                % AIC or not. If the new AIC is better, we change the best
                % parameters vector

                % We as well save the results in the all parameters matrix

                if G12_AIC < bestPars(6)
                    
                    allPars(m-1,:) = bestPars;
                    bestPars = [posParams, G12_pars, G12_AIC, parsVec(m)];

                else

                    allPars(m-1,:) = [posParams, G12_pars, G12_AIC, parsVec(m)];

                end
            
            end

        end

        % For posterior analysis of delay values
        best_tau2(bt2) = floor(4*bestPars(5))/4;
        bt2 = bt2+1;
        
        % For posterior analysis of best (delta, p) pairs
        if bestPars(7) == 2

            best_delta2(bd2) = bestPars(3);
            bd2 = bd2+1;

            best_p2(bp2) = bestPars(4);
            bp2 = bp2+1;

        end

        % Solution with optimal parameters, initialized with curve MV
        [tTG, xTG] = optimCurves(fitMean, TIV, ...
            [bestPars(3:4), floor(4*bestPars(5))/4], bestPars(1:2), ...
            [-3, day_NR+3+daysAfterTreat], hstep, day_NR-1, MVcurve);
    
        logxTG = log10(xTG(:,3)); % log10 viral load
        
        % We plot the best fitting in Figures 5 and 6
        figure(5)
        subplot(4, 4, mod(j-1,16)+1)

        plot(tTG, logxTG, '-', 'Color', colorsPars{bestPars(7)}, ...
            'LineWidth', 1.5)
        
        figure(6)
        subplot(4, 4, mod(j-1,16)+1)

        plot(tTG, logxTG, '-', 'Color', colorsPars{bestPars(7)}, ...
            'LineWidth', 1.5)

        for m = 1:2

            % Solution with non-optimal parameters, initialized with curve MV
            [tTG, xTG] = optimCurves(fitMean, TIV, ...
                [allPars(m,3:4), floor(4*allPars(m,5))/4], allPars(m,1:2), ...
                [-3, day_NR+3+daysAfterTreat], hstep, day_NR-1, MVcurve);
        
            logxTG = log10(xTG(:,3));
            
            % If AICs between best parameter tuple and this parameter tuple
            % are comparable (i.e., distance <= 2), then we also consider
            % the current parameter tuple as a best fitting
            if abs(allPars(m,6)-bestPars(6)) <= 2

                figure(5)
                subplot(4, 4, mod(j-1,16)+1)
        
                plot(tTG, logxTG, '-', 'Color', colorsPars{allPars(m,7)}, ...
                    'LineWidth', 1.5)

                best_tau2(bt2) = floor(4*allPars(m,5))/4;
                bt2 = bt2+1;

                if allPars(m,7) == 2

                    best_delta2(bd2) = allPars(m,3);
                    bd2 = bd2+1;

                    best_p2(bp2) = allPars(m,4);
                    bp2 = bp2+1;

                end

            end

            figure(6)
            subplot(4, 4, mod(j-1,16)+1)
    
            plot(tTG, logxTG, '-', 'Color', colorsPars{allPars(m,7)}, ...
                'LineWidth', 1.5)

        end

    end

end

% Create table 
filename = sprintf('TreatedGroupFittings.xlsx');
writetable(dataTableG12, filename);


%% Plots summarizing the results

% We recover the fitting parameters of Untreated Group
UGparams = readtable('UntreatedGroupFittings.xlsx');
UGbeta = table2array(UGparams(1:end-1, 2));
UGdelta = table2array(UGparams(1:end-1, 3));
UGp = table2array(UGparams(1:end-1, 4));
UGc = table2array(UGparams(1:end-1, 5));
UGR0 = table2array(UGparams(1:end-1, 6));
MVvalues = table2array(UGparams(end, 2:6)); % Corresponding to MV curve

xBP_UG = zeros(1,41); % Auxiliary zeroes vector to use for the box plots
% associated to Untreated Group
C = orderedcolors('gem');

% Box plots for {beta, delta, p, c, R0}
figure(7)

% beta
subplot(2, 3, 1)
boxchart(xBP_UG, UGbeta)
hold on
plot(xBP_UG, UGbeta, '.', 'MarkerSize', 10)
plot(0, MVvalues(1), 'diamond', 'MarkerFaceColor', 'black', ...
    'MarkerEdgeColor', 'black', 'MarkerSize', 7)
xlabel('\beta_{best}'), ylabel('ml/(copies x day)')

% delta
subplot(2, 3, 2)
boxchart(xBP_UG, UGdelta)
hold on
plot(xBP_UG, UGdelta, '.', 'MarkerSize', 10)
plot(0, MVvalues(2), 'diamond', 'MarkerFaceColor', 'black', ...
    'MarkerEdgeColor', 'black', 'MarkerSize', 7)
xlabel('\delta_{best}'), ylabel('1/day')

% p
subplot(2, 3, 3)
boxchart(xBP_UG, UGp)
hold on
plot(xBP_UG, UGp, '.', 'MarkerSize', 10)
plot(0, MVvalues(3), 'diamond', 'MarkerFaceColor', 'black', ...
    'MarkerEdgeColor', 'black', 'MarkerSize', 7)
xlabel('p_{best}'), ylabel('copies/(cells x day)')

% c
subplot(2, 3, 4)
boxchart(xBP_UG, UGc)
hold on
plot(xBP_UG, UGc, '.', 'MarkerSize', 10)
plot(0, MVvalues(4), 'diamond', 'MarkerFaceColor', 'black', ...
    'MarkerEdgeColor', 'black', 'MarkerSize', 7)
xlabel('c_{best}'), ylabel('1/day')

% R0
subplot(2, 3, 5)
boxchart(xBP_UG, UGR0, 'HandleVisibility', 'off')
hold on
plot(xBP_UG, UGR0, '.', 'MarkerSize', 10)
plot(0, MVvalues(5), 'diamond', 'MarkerFaceColor', 'black', ...
    'MarkerEdgeColor', 'black', 'MarkerSize', 7)
xlabel('R_0')
legend('Individual values', 'MV values')

sgtitle(['Box plots of the best parameters and corresponding R_0 for ' ...
    'Untreated Group'])


% Box plots for peak values and times, and clearance times
figure(8)

% Peak values
subplot(2, 2, 1)
boxchart(xBP_UG, peak_values(1:end-1))
hold on
plot(xBP_UG, peak_values(1:end-1), '.', 'MarkerSize', 10)
plot(0, peak_values(end), 'diamond', 'MarkerFaceColor', 'black', ...
    'MarkerEdgeColor', 'black', 'MarkerSize', 7)
xlabel('Peak value'), ylabel('log10 copies/ml')
title('A: Peak values')

% Peak times
subplot(2, 2, 2)
boxchart(xBP_UG, peak_times(1:end-1), 'HandleVisibility', 'off')
hold on
plot(xBP_UG, peak_times(1:end-1), '.', 'MarkerSize', 10)
plot(0, peak_times(end), 'diamond', 'MarkerFaceColor', 'black', ...
    'MarkerEdgeColor', 'black', 'MarkerSize', 7)
ylabel('days')

% Clearance times in the same subplot
boxchart(xBP_UG+1, clearance_times(1:end-1), 'BoxFaceColor', C(1,:), ...
    'MarkerColor', C(1,:), 'HandleVisibility', 'off')
hold on
plot(xBP_UG+1, clearance_times(1:end-1), '.', 'MarkerSize', 10, ...
    'MarkerEdgeColor', C(2,:), 'HandleVisibility', 'off')
plot(1, clearance_times(end), 'diamond', 'MarkerFaceColor', 'black', ...
    'MarkerEdgeColor', 'black', 'MarkerSize', 7, 'HandleVisibility', 'off')
xticks([0, 1])
xticklabels({'Peak time', 'Clearance time'})
legend('Individual values', 'MV values')
title('B: Peak and clearance times')

% Correlation between peak values and times
subplot(2, 2, 3)
plot(peak_times(1:end-1), peak_values(1:end-1), 'o', 'MarkerEdgeColor', ...
    C(2,:), 'MarkerFaceColor', C(2,:))
hold on
plot(peak_times(end), peak_values(end), 'diamond', 'MarkerEdgeColor', ...
    'black', 'MarkerFaceColor', 'black', 'MarkerSize', 7)
xlabel('Peak times (days)'), ylabel('Peak values (copies/ml)')

% Compute linear correlation
[pLinear, SLinear] = polyfit(peak_times(1:end-1), peak_values(1:end-1),1);
yLinear = polyval(pLinear, -1.5:8);
plot(-1.5:8, yLinear, '-', 'Color', C(3,:), 'LineWidth', 1.5)
R2Value = SLinear.rsquared;
sR2val = ['R2 value = ', num2str(R2Value)];
text(0, 6.4, sR2val)
title('C: Peak times and values correlation')

sgtitle('Peak values, peak times and clearance times for Untreated Group')


% Auxiliary zeroes vectors to use for the box plots associated to Treated
% Group
xd1 = zeros(1,length(best_delta1));
xd2 = zeros(1,length(best_delta2));
xp1 = zeros(1,length(best_p1));
xp2 = zeros(1,length(best_p2));

% Table to compare one MOA and two MOAs best delay (tau) values
tau1 = tabulate(best_tau1);
tau2 = tabulate(best_tau2);
tau2 = [0.25, 0, 0; tau2];


% Box plots for best delta and p values for one and two MOAs scenarios
figure(9)

% delta
subplot(2, 2, 1)
boxchart(xd1, best_delta1, 'BoxFaceColor', C(1,:), 'MarkerColor', C(1,:))
hold on
plot(xd1, best_delta1, '.', 'MarkerSize', 10, 'MarkerEdgeColor', C(2,:))
boxchart(xd2+1, best_delta2, 'BoxFaceColor', C(1,:), 'MarkerColor', C(1,:))
plot(xd2+1, best_delta2, '.', 'MarkerSize', 10, 'MarkerEdgeColor', C(2,:))
xticks([0 1])
xticklabels({'One MOA fits', 'Two MOAs fits'})
title('A: Box plots for best \delta_a values')
ylabel('1/day')
ylim([1, 3.5])

% p
subplot(2, 2, 2)
boxchart(xp1, best_p1, 'BoxFaceColor', C(1,:), 'MarkerColor', C(1,:))
hold on
plot(xp1, best_p1, '.', 'MarkerSize', 10, 'MarkerEdgeColor', C(2,:))
boxchart(xp2+1, best_p2, 'BoxFaceColor', C(1,:), 'MarkerColor', C(1,:))
plot(xp2+1, best_p2, '.', 'MarkerSize', 10, 'MarkerEdgeColor', C(2,:))
xticks([0 1])
xticklabels({'One MOA fits', 'Two MOAs fits'})
title('B: Box plots for best p_a values')
ylabel('copies/(cells x day)')
ylim([0, 0.8])

% Comparison delta vs. p for best (delta, p) tuples in the two MOAs 
% scenario
subplot(2, 2, 3)
plot(best_p2, best_delta2, 'o', 'MarkerEdgeColor', C(2,:), ...
    'MarkerFaceColor', C(2,:), 'HandleVisibility', 'off')
hold on
plot(fitMean(3), fitMean(2), 'diamond', 'MarkerFaceColor', 'black', ...
    'MarkerEdgeColor', 'black', 'MarkerSize', 7)
xlabel('p_a (copies/(cells x day))'), ylabel('\delta_a (1/day)')
legend('(p_{best}, \delta_{best})')
title('C: p_a vs. \delta_a')

sgtitle('\delta_a and p_a distributions for Treated Group')


% Frequencies of best delay values for one and two MOAs scenarios
figure(10)

bplot = bar([tau1(:,2).'; tau2(:,2).'], 'FaceColor', 'flat', ...
    'HandleVisibility', 'off');

for l = 1:7
    bplot(l).CData = [C(1,:); C(2,:)];
end

hold on
patch(NaN, NaN, C(1,:), 'DisplayName', 'One MOA fits')
patch(NaN, NaN, C(2,:), 'DisplayName', 'Two MOAs fits')
legend

xtick_coords = sort([bplot(1).XEndPoints, bplot(2).XEndPoints, ...
    bplot(3).XEndPoints, bplot(4).XEndPoints, bplot(5).XEndPoints, ...
    bplot(6).XEndPoints, bplot(7).XEndPoints]);
xtick_labels = {'0.25', '0.5', '0.75', '1', '1.25', '1.5', '1.75', ...
    '0.25', '0.5', '0.75', '1', '1.25', '1.5', '1.75'};
set(gca, 'XTick', xtick_coords, 'XTickLabel', xtick_labels);

title('Frequencies of best \tau values')
xlabel('days')
ylabel('Frequency')
ylim([0, 10])

%% Auxiliary functions

function err = treatment(UGparams, TIV, newParams, posParam, tSpan, ...
    hstep, dataValues, dataTime1, dataTime2, treatDay, UGcurves)
% Function to fit one parameter and the delay during treatment, computing
% the error between the fitted curve and patient data

% INPUTS:
% UGparams: basal parameters given by Untreated Group for initialization,
% (beta, delta, p, c)
% TIV: the TIV system
% newParams: the quantities we are fitting, (parameter, delay)
% posParam: the position of the parameter we are fitting, (beta=1, delta=2,
%  p=3, c=4)
% tSpan: time interval to compute solutions
% hstep: time step for integration
% dataValues: patient data
% dataTime1: times where patient data is measurable, > 1
% dataTime 2: times where patient data is unmeasurable, = 1
% treatDay: day that treatment starts
% UGcurves: the mean fitted curves for Untreated Group, (T,I,V)

% OUTPUTS:
% err: error between the log10 of the computed curve V and the patient data

    dataValues = dataValues(:);
    params = UGparams;

    for i = 1:length(posParam)
        
        % We substitute in the basal parameters the parameter we are
        % fitting
        params(posParam(i)) = newParams(i);

    end

    tau = newParams(end);
    
    % We need to obtain delay values that match the time values of our t
    % spans - we consider quarters of day, i.e., 6 hours
    tau = floor(4*tau)/4;
    
    % Define TIV model given the parameters
    TIVtreat = @(t,x) TIV(t,x,params);

    % Compute the solution since treatment starts effecting until 3 days
    % after patient has stopped receiving treatment
    [~, xTreat] = rungeKutta4(TIVtreat, [treatDay+tau, tSpan(2)], ...
        UGcurves((treatDay+tau-tSpan(1))/hstep + 1, :), hstep);

    % We divide the times depending on whether patient data is measurable 
    % or not
    tVec = tSpan(1):hstep:tSpan(2);
    dataTime1 = dataTime1(dataTime1 <= tSpan(2)); % Pick times inside our time span
    dataTime2 = dataTime2(dataTime2 <= tSpan(2)); % Pick times inside our time span
    dataValues = dataValues(1:length(dataTime1)); % Pick data inside our time span
 
    viralLoadTimes = ismember(tVec, dataTime1); % Measurable data
    noLoadTimes = ismember(tVec, dataTime2); % Unmeasurable data
    
    logV = log10([UGcurves(1:(treatDay+tau-tSpan(1))/hstep, 3); ...
        xTreat(:, 3)]); % log10 of viral load
    % Viral load is [MV curve until treatment starts effecting; new curve
    % during treatment effect]
    
    % Compute the error as the sum of:
    % 1. The norm of the difference between patient data and the computed
    % data
    % 2. Whenever there are times with unmeasurable data, a penalization to
    % any computed data above 1 in the days where patient data is
    % unmeasurable
    if ~isempty(tVec(noLoadTimes))
        
        data_nonMes = logV(noLoadTimes) - 1;
        data_nonMes(data_nonMes < 0) = zeros(length(data_nonMes(data_nonMes < 0)), 1);
        err = norm([logV(viralLoadTimes) - dataValues(dataTime1 <= tSpan(2)); ...
            data_nonMes]);

    else

        err = norm(logV(viralLoadTimes) - dataValues(dataTime1 <= tSpan(2)));
        
    end

end

function [tVec, x] = optimCurves(UGparams, TIV, newParams, posParam, ...
    tSpan, hstep, treatDay, UGcurves)
% Function to compute the fitted curve of invidiuals in Groups 1 and 2,
% given the fitted parameter and the delay

% INPUTS:
% UGparams: basal parameters given by Untreated Group for initialization,
% (beta, delta, p, c)
% TIV: the TIV system
% newParams: the quantities we are fitting, (parameter, delay)
% posParam: the position of the parameter we are fitting, (beta=1, delta=2,
% p=3, c=4)
% tSpan: time interval to compute solutions
% hstep: time step for integration
% treatDay: day that treatment starts
% UGcurves: the mean fitted curves for Untreated Group, (T,I,V)

% OUTPUTS:
% tVec: time vector of the solution
% x: curve V for treated patient

    params = UGparams;
    
    for i = 1:length(posParam)
        
        % We substitute in the basal parameters the parameter we are
        % fitting
        params(posParam(i)) = newParams(i);

    end

    tau = newParams(end);

    TIVtreat = @(t,x) TIV(t,x,params);
    [~, xTreat] = rungeKutta4(TIVtreat, [treatDay+tau, tSpan(2)],...
        UGcurves((treatDay+tau-tSpan(1))/hstep + 1, :), hstep);

    tVec = tSpan(1):hstep:tSpan(2);

    x = UGcurves(1:(treatDay+tau-tSpan(1))/hstep, :);
    x((treatDay+tau-tSpan(1))/hstep + 1:(tSpan(2)-tSpan(1))/hstep + 1,:) = xTreat;

end

function err = error_fun2(tspan, ti1, ti2, patient, x0, h, params)
% Function to compute the error when fitting the parameters of Untreated
% Group

% INPUTS:
% tspan: time interval to compute V
% ti1: times where data has been measured, > 1
% ti2: times where data was unmeasurable and set = 1
% patient: patient data at times ti1 and ti2
% x0: initial condition for curve V
% h: time step for integration
% params: parameters of TIV model

% OUTPUT:
% err: error between the log10 of the computed curve V and the data of the
% patient

    TIV = @(t,x) [- params(1)*x(1)*x(3);
        params(1)*x(1)*x(3) - params(2)*x(2);
        params(3)*x(2) - params(4)*x(3)];
    
    % Integration of the system via the Runge-Kutta 4 algorithm
    [t, x] = rungeKutta4(TIV, tspan, x0, h);

    % We extract the points in curve V such that we have patient data, and
    % depending on whether it is measurable or unmeasurable data
    x_data1 = x(ismember(t, ti1), 3).';

    if ~isempty(ti2)

        x_data2 = x(ismember(t, ti2), 3).';

    else

        x_data2 = 10;

    end

    % We compute the error as the sum of:
    % 1. The norm of the difference between patient data and the computed
    % data
    % 2. A penalization to any computed data above 1 in the days where
    % patient data is unmeasurable
    data_nonMes = log10(x_data2) - 1;
    data_nonMes(data_nonMes < 0) = zeros(1, length(data_nonMes(data_nonMes < 0)));
    err = norm([log10(x_data1) - patient, data_nonMes]);

end

function xdata = rungeKutta4aux(f, tspan, x0, h)
% Auxiliary function that returns only the curve of V

% INPUTS:
% f: ODE system to integrate
% tspan: time interval in which we integrate
% x0: initial condition of the system
% h: time step for integration

% OUTPUT:
% xdata: third solution curve
    
    % Integration of the system via the Runge-Kutta 4 algorithm
    [~, x] = rungeKutta4(f, tspan, x0, h);
    
    % Third column corresponds to V
    xdata = x(:, 3);

end

function [t, x] = rungeKutta4(f, tspan, x0, h)
% Runge-Kutta 4 integration scheme

% INPUTS:
% f: ODE system to integrate
% tspan: time interval in which we integrate
% x0: initial condition of the system
% h: time step for integration

% OUTPUTS:
% t: time span of the solution
% x: numerical solution

    t = tspan(1):h:tspan(2);
    n = length(t);
    x = zeros(n, length(x0));
    x(1, :) = x0;
    
    for i = 1:n-1
        k1 = h * f(t(i), x(i, :));
        k2 = h * f(t(i) + h/2, x(i, :) + k1(:).'/2);
        k3 = h * f(t(i) + h/2, x(i, :) + k2(:).'/2);
        k4 = h * f(t(i) + h, x(i, :) + k3(:).');
        x(i+1, :) = x(i, :) + (k1(:) + 2*k2(:) + 2*k3(:) + k4(:)).' / 6;
    end

end