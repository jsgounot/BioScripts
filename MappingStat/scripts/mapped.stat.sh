# @Author: jsgounot
# @Date:   2020-12-09 11:18:39
# @Last Modified by:   jsgounot
# @Last Modified time: 2021-03-10 14:52:29

mapped=$(samtools view -c -F 4 $1)
unmapped=$(samtools view -c -f 4 $1)
sumvalue=$(($mapped+$unmapped))
prc=$(bc <<< "scale=4;$mapped/$sumvalue")
echo "mapped:$mapped"
echo "unmapped:$unmapped"
echo "sum:$sumvalue"
echo "prc:$prc"
