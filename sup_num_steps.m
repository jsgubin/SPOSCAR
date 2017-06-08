clear all;
global  loop max_iterations
loop=1;
max_iterations=10000;
data_flag=3;
d2=[0.5 1 2];
sizet=[ 5000  10000 15000 20000];
sizet=[ 1500  3000 4500 6000];
sizet=[ 800  1600 2400 3200];
length1=length(d2);
length2=length(sizet);
time_FastOSCAR=cell(length1,length2);
time_Path=cell(length1,length2);
Iterations=cell(length1,length2);
CallFastOscar=cell(length1,length2);
Singularities=cell(length1,length2);
load tmp11;
for i=3:length1 % d2
    local_d2 = d2(i);
    d=[1 local_d2];
    for j=2:length2                 % size_training
        size_training=sizet(1,j);
        [tmp1,tmp2,tmp3,tmp4,tmp5]=num_steps(data_flag,size_training,d);
        time_FastOSCAR{i,j}=tmp1;
        time_Path{i,j}=tmp2;
        Iterations{i,j}=tmp3;
        CallFastOscar{i,j}=tmp4;
        Singularities{i,j}=tmp5;
        save tmp11;
    end
end
save  lv_results;