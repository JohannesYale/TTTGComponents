% Copyright (c) 2024 Johannes Sieberer
% 
% Permission is hereby granted, free of charge, to any person obtaining a copy
% of this software and associated documentation files (the "Software"), to deal
% in the Software without restriction, including without limitation the rights
% to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
% copies of the Software, and to permit persons to whom the Software is
% furnished to do so, subject to the following conditions:
% 
% The above copyright notice and this permission notice shall be included in all
% copies or substantial portions of the Software.
% 
% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
% IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
% FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
% AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
% LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
% OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
% SOFTWARE.

%% Get patient files, intial results
listing = dir ("Dataset/*.csv");
cases = struct2table(listing).name;
results.TTTG = zeros(height(cases),1);
results = struct2table(results);
m = 0;
%% Iterate through dataset files
for n = 1:66
casen = cases(n-floor((m+1)/2));
data = importfile(strcat("Dataset/",casen)); %ImportData
    if(height(data) <= 23) %For unilateral knees
        %% Get Points (x,y,y)
        tibLatPostCon =  table2array(data(17,10:12)); %tibia lateral posterior condyle
        tibMedPostCon = table2array(data(18,10:12)) ; %tibia medial posterior condyle
        tibMedIterCondTub = table2array(data(16,10:12));  %tibia medial intercondylar tubercle
        tibLatIterCondTub = table2array(data(15,10:12));  %tibia lateral intercondylar tubercle
        tibTuber = table2array(data(20,10:12)); %Tibial tuberosity
        tibShaftCen =  table2array(data(19,10:12)); %Tibia shaft center
        femPostCondLat = table2array(data(5,10:12)); %femur posterior lateral condyle
        femPostCondMed = table2array(data(6,10:12)); %femur posterior medial condyle
        troch = table2array(data(22,10:12)); %trochlea groove
        results.Cases(n) = casen;
    else 
        if(mod(m,2) == 0) %Left knee bilateral
        tibLatPostCon =  table2array(data(33,10:12));
        tibMedPostCon = table2array(data(34,10:12)) ;
        tibMedIterCondTub = table2array(data(31,10:12));
        tibLatIterCondTub = table2array(data(29,10:12));
        tibTuber = table2array(data(39,10:12));
        tibShaftCen =  table2array(data(37,10:12));
        femPostCondLat = table2array(data(9,10:12));
        femPostCondMed = table2array(data(10,10:12));
        troch = table2array(data(43,10:12));
        m = m+1;
        results.Cases(n) = strcat(casen,"L");
    else %Right knee bilateral
        tibLatPostCon =  table2array(data(35,10:12));
        tibMedPostCon = table2array(data(36,10:12)) ;
        tibMedIterCondTub = table2array(data(30,10:12));
        tibLatIterCondTub = table2array(data(32,10:12));
        tibTuber = table2array(data(40,10:12));
        tibShaftCen =  table2array(data(38,10:12));
        femPostCondLat = table2array(data(11,10:12));
        femPostCondMed = table2array(data(12,10:12)); 
        troch = table2array(data(44,10:12));
        results.Cases(n) = strcat(casen,"R");
        m=m+1;
        end
    end
    
    %% Calculate TTTG, translational TTTG and Rotational TTTG

    % Get line tib shaftcenter to intercondylar midpoint
    tibMidPla = (tibMedIterCondTub + tibLatIterCondTub)/2;
    tibShaftMid = tibMidPla - tibShaftCen; %Tibial longitudinal axis

    %Get lines femur tibial rot
    femPostCondLine = femPostCondLat -femPostCondMed;
    femCondLineProj = femPostCondLine - dot(femPostCondLine,tibShaftMid)*tibShaftMid/norm(tibShaftMid)^2; 
    femCondLineProjNorm  = femCondLineProj/norm(femCondLineProj);%Projected femoral condyle line normalized
    tibPostCondLine = tibLatPostCon - tibMedPostCon;
    tibCondLineProj = tibPostCondLine - dot(tibPostCondLine,tibShaftMid)*tibShaftMid/norm(tibShaftMid)^2; 
    tibCondLineProjNorm =  tibCondLineProj/norm(tibCondLineProj); %Projected tibial condyle line normalized
    
    % Tibiofemoral rotation 
    rotVector = cross(femCondLineProj,tibCondLineProj); %Assume knee is generally superior inferior aligned (third component pointing up and down)
    results.ROTTF(n) = signedAngleTwo3DVectors(femCondLineProj,tibCondLineProj,rotVector(:) .* sign(rotVector(3)).* sign(femCondLineProj(1)),1)*180 / pi;

    %3D TT-TG tt to trochlea2 midpoint
    tttg = tibTuber - troch;
    results.TTTG(n) = dot(femCondLineProjNorm, tttg);

    %Tuberosity distance
    tttp =  tibMidPla - tibTuber; %tibial plateau to tibial tuberosity
    tubDist = -dot(tibCondLineProjNorm, tttp); %Tuberosity distance

    %Trochlea distance
    tgtp =  troch - tibMidPla; %tibial plateau to tibial tuberosity
    tgdist = -dot(femCondLineProjNorm, tgtp);  %Trochlea groove distance

    %Translational and rotational TTTG
    results.TransTTTG(n) = tgdist + tubDist;
    results.RotDist(n) = results.TTTG(n) -  results.TransTTTG(n);
end
