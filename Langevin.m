function [L_r2Mag,dL_r2Mag] = Langevin(r2_Mag)

L_r2Mag=coth(r2_Mag)-1./r2_Mag; 
dL_r2Mag=1./(r2_Mag.^2)-(coth(r2_Mag)).^2+1;

    
end
