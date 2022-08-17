function [P, P2, sontemp1] = BayesEstimatorL(dist, Pu, w1, w2, sigma2) 
    
    if w2<1
        [Nei,D] = size(dist);    
        mon2 = ((1-w1)/w1) * ((2*pi*sigma2)^(Nei*D/2)*Pu^Nei) ;
        sontemp = dist / (2*sigma2);
        sontemp1 = sum(sontemp,2); 
        sonexp1 = w2*exp(-sontemp1);      sonexp2 = (1-w2)*(2*pi*sigma2^(D/2))*Pu;  % ¸ÅÂÊÏòÁ¿
        P2 = sonexp1./(sonexp1+sonexp2);
        son = prod(sonexp1+sonexp2);
        P = son/(son+mon2);

    else  % 
        Nei = w2;
        mon2 = ((1-w1)/w1) * ((2*pi*sigma2)^(Nei*D/2)*Pu^Nei) ;
        sontemp = dist / (2*sigma2);
        sontemp1 = sort(sum(sontemp,2)); 
        sonexp1 = exp(-sontemp1(1:Nei)); 
        son = prod(sonexp1);
        P = son/(son+mon2);     P2 = zeros(Nei,1);
    end