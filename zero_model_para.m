function [model_disch] = zero_model_para(R,PET,lag)
% R = rainfall time series  (mm/day)
% PET = potential evapotranspiration  time series in mm/day
% lag =  use the value between 0 to 2.  ( to consider channel travel time, see original paper for details) 

a=0.02 ;   
%% all parameter values
            W=zeros(1,1);
            H=zeros(1,1);
            
            len = size(R,1);
            % Calculate W
            for iii= 1:len-1
                if R(iii,1)> PET(iii,1)
                    W(iii,1) = R(iii,1)-PET(iii,1);
                else
                    W(iii,1) = 0;
                end
            end
            % Calculate H
            for i= 1:len-1
                if PET(i,1)> R(i,1)
                    H(i,1) = PET(i,1)-R(i,1);
                else
                    H(i,1) = 0;
                end
            end
            
            % Calculate FW
            FW = zeros(len,1);
            for i = 1:len-1
                if i<366
                    FW(i,1)=  nan;
                else
                    count = 0;
                    for j = 1:365
                        FW(i,1)= FW(i,1)+ W(i-j+1,1)*exp(-count*a);
                        count = count + 1;
                    end
                end
            end
  %---------------------------**Important**----------------------------------------------
   % in some cases FW is coming as zero. so the 'phi' value which is ratio
   % of 'FH' by 'FW'  comes as infinity. so zero of FW by 0.00001 to avoid
   % the nan value. 
     for nn=1:size(FW,1)
         if FW(nn,1)==0
            FW(nn,1)=0.0000001 ;
         end
     end
   
   %-----------------------------END of-----------------------------------
            % Calculate FH
            FH = zeros(len,1);
            for i = 1:len-1
                if i<366
                    FH(i,1)=  nan;
                else
                    count = 0;
                    for j = 1:365
                        FH(i,1)= FH(i,1)+ H(i-j+1,1)*exp(-count*a);
                        count = count + 1;
                    end
                end
            end
            % Calculate FW/FW+FH
%             FW_ratio = FW./(FW+FH);
%             FH_ratio = FH./(FW+FH);
            
            % Calculate Discharge
            
            model_disch = zeros(len,1);
            % Calculate
            phi = FH./FW;
            IR = zeros(length(phi),1);
            for i = 366:length(phi)-1
                IR(i,1) = W(i,1)* (1-((phi(i,1)*tanh(phi(i,1)^(-1)))*(1-exp(-phi(i,1))))^0.5);
            end
            for i=365:len-2
                if i<=800
                    
                    for k=365:i
                        
                        model_disch(i+lag,1)= model_disch(i+lag,1) + IR(k,1)*0.4/((i-k)*0.4+1)/((i-k+1)*0.4+1);
                    end
                else
                    for k=i-365:i
                        
                        model_disch(i+lag,1)= model_disch(i+lag,1) + IR(k,1)*0.4/((i-k)*0.4+1)/((i-k+1)*0.4+1);
                    end
                end
            end
            model_disch(1:366,1) = NaN ;
end

