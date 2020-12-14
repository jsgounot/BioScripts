# @Author: jsgounot
# @Date:   2020-12-09 11:18:39
# @Last Modified by:   jsgounot
# @Last Modified time: 2020-12-09 11:31:55

input=$(realpath $1)
mapped=$(samtools view -c -F 4 $input)
unmapped=$(samtools view -c -f 4 $input)
sumvalue=$(($mapped+$unmapped))
prc=$(bc <<< "scale=4;$mapped/$sumvalue")
echo "mapped:$mapped"
echo "unmapped:$unmapped"
echo "sum:$sumvalue"
echo "prc:$prc"