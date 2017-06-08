function [out]=main(data_flag,size,dir)
    global KTYPE KSCALE fake_zero size_training fake_zero_linprog NumSingular  time_path
    % KTYPE = 6;
    fake_zero=10^-8;
    fake_zero2=10^-8;
    NumSingular=0;
    time_batch =0;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [x,y]=read_data(data_flag);
    [training_x,training_y,test_x,test_y]=random_select(x,y,size);
    save aa;
    load aa;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    uper = SPOSCAR.ComputUpper(training_x,training_y,dir);
    range_zeta = [0.1 uper];
    a1 = tic ;
    out=SPOSCAR(training_x,training_y,dir,range_zeta);
    time_path=toc(a1);
%     a=OSCAR(training_x(:,[1:400]),training_y,dir);
%     a=OSCAR(training_x,training_y,dir);
end
