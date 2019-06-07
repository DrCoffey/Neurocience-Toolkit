function [ druglevel, time] = pharmacokineticsV3(varargin)
%% This fucntion calculates cocaine drug level (Molar Mass of Cocaine - 339.81g/MOLE - SigmaAldrich)
% Written By Kevin Coffey (Rutgers) & Olivia Kim (UPenn)
% INPUTS: infusions(ms_inf_onset,ms_inf_offset)
%         weight(body weight in Kg);
%         type();
%             (1) = IV - First Order (Eliminatoin Only; Root, 2011)
%             (2) = IV - Estimated Brain Level (2 compartmaent Model; Pan 1991; Roberts, 2013)
%             (3) = IV - Estimated Blood Level (2 compartmaent Model; Pan 1991)
%             (4) = IP - Estimated Brain Level (2 compartmaent Model; Pan 1991)
%             (5) = IP - Estimated Blood Level (2 compartmaent Model; Pan 1991)
%
%         sesmins(Length of Session in Minutes)
%         makeFig(1 make figure / 0 no figure)
%
% OUTPUT: Drug Level Duh!

p = inputParser;
val = @(x)validateattributes(x,{'numeric'},{'scalar'});
p.addRequired('infusions', @(x)validateattributes(x,{'double'},{'ncols',2})); % [infusionOn infusionOff] (ms)
p.addParameter('type', 2, @(x)mustBeMember(x,[1,2,3,4,5])); % Equation type
p.addParameter('duration', [], val); % Duration of the session in minutes

% Either use weight, mg/mL, and mg/Kg, OR use mg/Kg/sec
p.addParameter('weight', 500, val); % Rat weight in grams
p.addParameter('mg_mL', [], val); % Cocaine concentration (mg/mL)
p.addParameter('mL_S', [], val); % pump flow rate (mL/sec)
p.addParameter('mg_Kg_s', 0.366 / 3.156, val); % Cocaine mg/Kg/sec
        
p.parse(varargin{:});
%% Calculate uM per Kg per second
molarMass = .33981; % Cocaine HCL (339.81g/MOLE or .33981mg/uMOLE)
% If mg/mL, mL/sec, and weight are given, use this to calculate mg/Kg/second
if length([p.Results.mg_mL, p.Results.mg_mL, p.Results.weight]) == 3
    MGKGDOSE = p.Results.mg_mL * p.Results.mL_S; % mg/sec = mg/mL * mL/sec
    MGKGDOSE = 1000 * MGKGDOSE / p.Results.weight; % mg/Kg/sec = mg/sec / (mg * 1000)
else % use the value supplied
    MGKGDOSE = p.Results.mg_Kg_s; % Cocaine mg/Kg/sec
end
uMDOSE = MGKGDOSE/molarMass; % uM per Kg per second 

%% Length of the session
SESSIONLENGTH = p.Results.duration;
if isempty(SESSIONLENGTH)
    SESSIONLENGTH = max(max(p.Results.infusions)) / 60000;
end
SESSIONLENGTH = round(SESSIONLENGTH);

% Infusion start, end, and duration
infusionStart = p.Results.infusions(:,1) / 1000; % Infusion Start in 1s Units
infusionEnd   = p.Results.infusions(:,2) / 1000; % Infusion End in 1s Units
infusionDuration = infusionEnd - infusionStart; % Infusion Duration in 1s Units

        
switch p.Results.type
    %% %%%%%%%%%%%%%%% TYPE 1: First Order Calculation %%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%  Output Concentration in mg/kg  %%%%%%%%%%%%%%%%%%%%%%%%
    case 1
        % Variables For the (Type 1) First Order Pharmacokinetics Equation
        HALFLIFE = 24; % Enther Half Life in Minutes
        K = log(2)/HALFLIFE; % Calculate Elimination Paramenter
        infusioncheck=zeros(SESSIONLENGTH*60,2); % Logical Array Containing Infusion Logic
        timespan=(0:(SESSIONLENGTH*60)-1)'; % Time Array
        p=[];
        for i=1:length(infusionStart) % Calc Infusion Logical Array (Infusion or Not)
            if infusionStart(i)>=sesmins*60;
                continue %Ignore infusions after session end
            else
                % Calc Infusion Logical Array (Infusion or Not)
                p(end+1,1)=find((timespan(:,1)==floor(infusionStart(i))));
            end
        end
        infusioncheck(p)=1;
        for i=1:length(infusionStart)
            % Calc Infusion Dur (Infusion or Not)
            infusioncheck(floor(infusionStart(i))+1,2)=infusionDuration(i);
        end
        druglevel = 0;
        for i=1:length(infusioncheck)-1
            druglevel(i+1)=(druglevel(i,1)+(MGKGDOSE*infusioncheck(i+1,2))*infusioncheck(i+1))*exp(-K/(60));
        end
     
        %% %%%%%%%%%%%%%%% TYPE 2: Intravenus Estimated Brain Level %%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%  Output Concentration in uMole/L or uM   %%%%%%%%%%%%%%
    case 2
        % Variables For the (Type 2) 2 Compartment - Estimated Brain Levels
        k12 = 0.233; % Pan & Justice 1990
        k21 = 0.212; % Pan & Justice 1990
        kel = 0.294; % Pan & Justice 1990
        ALPHA = 0.5*((k12 + k21 + kel)+sqrt((k12 + k21 + kel)^2-(4*k21*kel))); %%%% Represent the redistribution of Cocaine. Calculated using Pan et al. 1991 eqn. 0.06667 with Lau and Sun 2002 values for 0.5 mg/kg dose (see paper for justification)
        BETA = 0.5*((k12 + k21 + kel)-sqrt((k12 + k21 + kel)^2-(4*k21*kel))); %%%% Represent the Elimination of Cocaine. Calculated using Pan et al. 1991 eqn. 0.0193 here with Lau and Sun 2002 values for 0.5 mg/kg dose
        VOLUME = .15; %%%% Brain Apperant Volume of distribution in L per kg
        % K_FLOW = .233; %%%% Represents the flow between the two compartments
        inf_dl=zeros(length(infusionDuration),SESSIONLENGTH*60); % Pre allocate array
        % Calculate Drug Level for Each Infusion seperately
        for i = 1:length(infusionDuration)
            for j = round(infusionStart(i,1))+1:(SESSIONLENGTH*60)
                inf_dl(i,j)=(((uMDOSE*infusionDuration(i,1))*(k12))/(VOLUME*(ALPHA-BETA)))*(exp(-BETA*((j-round(infusionStart(i,1)))/60))-exp(-ALPHA*((j-round(infusionStart(i,1)))/60)));
            end
        end
        druglevel=sum(inf_dl); % Sum Individual Infusion Drug Levels
       
        %% %%%%%%%%%%%%%%% TYPE 3: Estimated Blood Level %%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%     Concentration in uMole/L or uM       %%%%%%%%%%%%%%
    case 3
        % Variables For the (Type 3) - Estimated Blood Levels Equation
        k12 = 0.233; % Pan & Justice 1990
        k21 = 0.212; % Pan & Justice 1990
        kel = 0.294; % Pan & Justice 1990
        ALPHA = 0.5*((k12 + k21 + kel)+sqrt((k12 + k21 + kel)^2-(4*k21*kel))); %%%% Represent the redistribution of Cocaine. Calculated using Pan et al. 1991 eqn. 0.06667 with Lau and Sun 2002 values for 0.5 mg/kg dose (see paper for justification)
        BETA = 0.5*((k12 + k21 + kel)-sqrt((k12 + k21 + kel)^2-(4*k21*kel))); %%%% Represent the Elimination of Cocaine. Calculated using Pan et al. 1991 eqn. 0.0193 here with Lau and Sun 2002 values for 0.5 mg/kg dose
        VOLUME = .120; %%%% Apperant Volume of blood (Volume of Distrobution in L/kg). From Lau & Sun 2002
        % K_FLOW = .233; %%%% Represents the flow between the two compartments
        inf_dl=zeros(length(infusionDuration),SESSIONLENGTH*60); % Pre allocate array
        % Calculate Drug Level for Each Infusion seperately
        for i = 1:length(infusionDuration)
            for j = round(infusionStart(i,1))+1:(SESSIONLENGTH*60)
                inf_dl(i,j)=(((uMDOSE*infusionDuration(i,1)))/(VOLUME*(ALPHA-BETA)))*(((k12-BETA)*exp(-BETA*((j-round(infusionStart(i,1)))/60)))-((k12-ALPHA)*exp(-ALPHA*((j-round(infusionStart(i,1)))/60))));
            end
        end
        druglevel=sum(inf_dl); % Sum Individual Infusion Drug Levels
        
        
        %% %%%%%%%%%%%%%%% TYPE 4: Intraparitoneal Estimated Brain Level %%%%%%%%%
        %%%%%%%%%%%%%%%%%%     Concentration in uMole/L or uM       %%%%%%%%%%%%%%
    case 4
        % Variables For the (Type 4) 2 Compartment - IP Estimated Brain Levels
        k12 = 0.233; % Pan & Justice 1990
        k21 = 0.212; % Pan & Justice 1990
        kel = 0.294; % Pan & Justice 1990
        ALPHA = 0.5*((k12 + k21 + kel)+sqrt((k12 + k21 + kel)^2-(4*k21*kel))); %%%% Represent the redistribution of Cocaine. Calculated using Pan et al. 1991 eqn. 0.06667 with Lau and Sun 2002 values for 0.5 mg/kg dose (see paper for justification)
        BETA = 0.5*((k12 + k21 + kel)-sqrt((k12 + k21 + kel)^2-(4*k21*kel))); %%%% Represent the Elimination of Cocaine. Calculated using Pan et al. 1991 eqn. 0.0193 here with Lau and Sun 2002 values for 0.5 mg/kg dose
        VOLUME = .15; %%%% Brain Apperant Volume of distribution in L per kg
        F=.8581; % Derived from 30mg/kg estimated parameter (F*88.28uMol/kg)/(Volime Distribution Blood l/kg)-(Pan & Justice 1990)
        kA= .0248; % Pan & Justice 1990
        % K_FLOW = .233; %%%% Represents the flow between the two compartments
        inf_dl=zeros(length(infusionStart),SESSIONLENGTH*60);% Pre allocate array
        % Calculate Drug Level for Each Infusion seperately
        for i = 1:length(infusionStart)
            for j = round(infusionStart(i,1))+1:(SESSIONLENGTH*60)
                inf_dl(i,j)=(F*uMDOSE*kA*k21/VOLUME)*...
                    ( (exp(-ALPHA*((j-round(infusionStart(i,1)))/60)))/((kA-ALPHA)*(BETA-ALPHA)) ...
                    + (exp(-BETA*((j-round(infusionStart(i,1)))/60)))/((kA-BETA)*(ALPHA-BETA)) ...
                    + (exp(-kA*((j-round(infusionStart(i,1)))/60)))/((ALPHA-kA)*(BETA-kA)) );
            end
        end
        if length(inf_dl(:,1))>1
            druglevel=sum(inf_dl); % Sum Individual Infusion Drug Levels
        else
            druglevel=inf_dl;
        end
        
        %% %%%%%%%%%%%%%%% TYPE 5: Intraparitoneal Estimated Blood Level %%%%%%%%
        %%%%%%%%%%%%%%%%%%     Concentration in uMole/L or uM       %%%%%%%%%%%%%
    case 5
        % Variables For the (Type 5) 2 Compartment - IP Estimated Blood Levels
        k12 = 0.233; % Pan & Justice 1990
        k21 = 0.212; % Pan & Justice 1990
        kel = 0.294; % Pan & Justice 1990
        ALPHA = 0.5*((k12 + k21 + kel)+sqrt((k12 + k21 + kel)^2-(4*k21*kel))); %%%% Represent the redistribution of Cocaine. Calculated using Pan et al. 1991 eqn. 0.06667 with Lau and Sun 2002 values for 0.5 mg/kg dose (see paper for justification)
        BETA = 0.5*((k12 + k21 + kel)-sqrt((k12 + k21 + kel)^2-(4*k21*kel))); %%%% Represent the Elimination of Cocaine. Calculated using Pan et al. 1991 eqn. 0.0193 here with Lau and Sun 2002 values for 0.5 mg/kg dose
        VOLUME = .12; %%%% Apperant Volume of blood (Volume of Distrobution in L/kg). From Lau & Sun 2002
        F=.3915; % Derived from 30mg/kg estimated parameter (F*88.28uMol/kg)/(Volime Distribution Blood l/kg)-(Pan & Justice 1990)
        kA= .0248; % Pan & Justice 1990
        % K_FLOW = .233; %%%% Represents the flow between the two compartments
        inf_dl=zeros(length(infusionStart),SESSIONLENGTH*60); % Pre allocate array
        % Calculate Drug Level for Each Infusion seperately
        for i = 1:length(infusionStart)
            for j = round(infusionStart(i,1))+1:(SESSIONLENGTH*60)
                inf_dl(i,j)=(F*uMDOSE*kA/VOLUME)*...
                    ( (((k21-ALPHA)/((kA-ALPHA)*(BETA-ALPHA)))*exp(-ALPHA*((j-round(infusionStart(i,1)))/60))) ...
                    + (((k21-BETA)/((kA-BETA)*(ALPHA-BETA)))*exp(-BETA*((j-round(infusionStart(i,1)))/60))) ...
                    + (((k21-kA)/((ALPHA-kA)*(BETA-kA)))*exp(-kA*((j-round(infusionStart(i,1)))/60))) );
            end
        end
        if length(inf_dl(:,1))>1
            druglevel=sum(inf_dl); % Sum Individual Infusion Drug Levels
        else
            druglevel=inf_dl;
        end
        
end
time = {[1:length(druglevel)]./60};
druglevel = {druglevel};

