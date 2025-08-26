#!/bin/bash

rm -f subjob.sh
ii=1001
while [ $ii -le 6000 ] #total 
do

cp -f test.py test-${ii}-1.py
sed -i "148,159s/seed/seed${ii}/g" test-${ii}-1.py
sed -i "160,200s/seed/seed${ii}-1/g" test-${ii}-1.py
sed -i "143s/ch2oo-2/ch2oo-2-${ii}/g" test-${ii}-1.py

nr=$((10 * ( ii - 1 ) + 1))
cp -f test.str test-seed${ii}.str
sed -i s:%seed%:${nr}:g test-seed${ii}.str

cp -f submit.sh sub-${ii}.sh
sed -i s:%seed%:${ii}:g sub-${ii}.sh

echo "sbatch sub-${ii}.sh" >> subjob.sh
for jj in $(seq 2 10);do
   cp test-a.py test-${ii}-${jj}.py
   kk=`expr $jj - 1`
   sed -i "143s/ch2oo-2/ch2oo-2-${ii}/g" test-${ii}-${jj}.py
   sed -i "159s/seed/seed${ii}/g" test-${ii}-${jj}.py
   sed -i "160,200s/seed/seed${ii}-${jj}/g" test-${ii}-${jj}.py
done
echo $ii

ii=`expr $ii + 1`
done

chmod +x subjob.sh

