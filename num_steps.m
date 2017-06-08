function [time_FastOSCAR,time_Path,Iterations,CallFastOscar,Singularities]=num_steps(data_flag,size_training,d)
    global loop time_batch time_path NumSingular
    time_FastOSCAR=[];
    time_Path=[];
    Iterations=[];
    CallFastOscar=[];
    Singularities=[];
    for i=1:loop
        flag=1;
        while flag==1
            try
                [out]=main(data_flag,size_training,d);
                flag=0;
            catch exception
                flag=1;
            end
        end
        time_FastOSCAR=[time_FastOSCAR;time_batch];
        time_Path=[time_Path;time_path];
        Iterations=[Iterations;out.Steps];
        CallFastOscar=[CallFastOscar;out.NumSubProblem];
        Singularities=[Singularities;NumSingular];
    end
end