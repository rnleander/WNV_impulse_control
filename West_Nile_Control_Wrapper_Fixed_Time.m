function [tau,ul,ua,X,J,J_comp,final_treatment_time,X0,T] = West_Nile_Control_Wrapper_Fixed_Time(N,Tf,Obj_type,larvicide)
%This code optimizes the control level at fixed application times
%in order to minimize a user selected objective functional. 
%Obj_type sets the objective functional. Vector control: 1, Disease
%control:2, Host preservation: 3. 
%larvicide_type: 1=long-lasting s-methorpene briquet, 2=VectoBac

if Obj_type==1
    [tau,ul,ua,X,J,J_comp,final_treatment_time,X0,T] = West_Nile_Control_Fixed_Times(N,Tf,larvicide);
end
if Obj_type==2
    [tau,ul,ua,X,J,J_comp,final_treatment_time,X0,T] = West_Nile_Control_Disease_Control_Fixed_Time(N,Tf,larvicide);
end
if Obj_type==3
    [tau,ul,ua,X,J,J_comp,final_treatment_time,X0,T] = West_Nile_Control_Host_Preservation_Fixed_Time(N,Tf,larvicide);
end


end
