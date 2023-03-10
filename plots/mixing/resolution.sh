flavor=('e' 'mu' 'tau')
operator=('U')

for i in ${operator[*]}
do
	for j in ${flavor[*]}
	do
		pdf2ps -r1400 "$i""$j"N_majorana.pdf
		mv "$i""$j"N_majorana.ps low_res/"$i""$j"N_majorana.ps
		ps2pdf low_res/"$i""$j"N_majorana.ps low_res/"$i""$j"N_majorana.pdf
		rm -r low_res/"$i""$j"N_majorana.ps
	done
done