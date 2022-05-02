img1 = imread('N:\citymapper_70\images\029_065_id3188c81615_124544_Backward.tif');
img2 = imread('N:\citymapper_70\images\029_056_id3179c82579_124523_Nadir.tif');

copts = load('N:\citymapper_70\res\Point_clouds\temp_folder_cluster\029_065_id3188c81615_124544_Backward\1.txt');

len = length(copts(:,1));

figure;
imshow(img2);
hold on;
plot(copts(:,1),copts(:,2),'*');
hold on;
text(copts(:,1),copts(:,2),num2str([1:len]'));


figure
imshow(img1);
hold on;
plot(copts(:,3),copts(:,4),'*');
hold on;
text(copts(:,3),copts(:,4),num2str([1:len]'));

