path=/projects/micb405/project1/Team9/project_2
cd $path/ffn

mags=(314 73 16 69 210 289 292 294 45 287 168 14 148 157 175 229 250 96)

for mag in ${mags[@]}
do
  cp $path/prokka/${mag}/*_bacteria.ffn $path/ffn/
done

for f in $path/ffn/*.ffn
do
  cat $f >> total.ffn
done
