% Oscar2DPath (Two-dimensional Solution Path Algorithm for Sparse Regression with Automatic Feature Grouping) demo 
% 
% Bin Gu, Nanjing Univerity of Information Secience and Technolgy
% 9/4/2016
% 
clear all;
data_flag=3;  % 1-4
size_training=3300; % size of dataset by randomly selecting
direction = [1 0.5];
[out]=main(data_flag,size_training,direction);  % 
save lv_results;
