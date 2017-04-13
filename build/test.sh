make
./ExtendedKF ../data/sample-laser-radar-measurement-data-1.txt output_1.txt > debug.txt
#python visualize.py -f output_1.txt &

./ExtendedKF ../data/sample-laser-radar-measurement-data-2.txt output_2.txt > debug2.txt
#python visualize.py -f output_2.txt &
echo "Data set 1 "
grep RMSE -A 5 debug.txt
echo "Data set 2 "
grep RMSE -A 5 debug2.txt
