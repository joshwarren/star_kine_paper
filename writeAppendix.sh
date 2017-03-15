#!/bin/bash
echo '\section{Zero moment Maps}'
echo '    \label{sec:0maps}'
for c in stellar OIII Hbeta Hgamma Hdelta NI 
do
	echo '    \subsection{'${c^}' maps}'
	echo '        \label{subsec:'$c'maps}'
	


	if [ $c = stellar ]
	then 
		comp='img'
	else
		comp='img eqW'
	fi

	if [ $c = OIII ]
	then
		d='[OIII]5007d'
		e='[OIII] '
	elif [ $c = NI ]
	then
		d='[NI]d'
		e='[NI] '
	else
		d=$c
		if [ $c = Hbeta ]
		then
			e='H$_\mathrm{\beta}$ '
		elif [ $c = Hgamma ]
		then
			e='H$_\mathrm{\gamma}$ '
		elif [ $c = delta ]
		then
			e='H$_\mathrm{\delta}$ '
		fi
	fi

	for p in $( echo $comp)
	do
		echo '        \begin{figure*}'
		echo '            \centering'
		# alphabetical 
		# for g in eso443-g024 ic1459 ic1531 ic4296 ngc0612 ngc1399 ngc3100 ngc3557 ngc7075 pks0718-34
		# rotational velocity
		for g in ngc0612 ngc3557 ngc3100 ic1459 pks0718-34 ic4296 ngc7075 ic1531 ngc1399 eso443-g024     
		do
			echo '            \includegraphics[width=0.245\textwidth]{plots/'$g'_'$d'_'$p'.png}'
		done

		if [ $p = img ]
		then 
			cap='Image'
		elif [ $p = eqW ]
		then 
			cap='Equivelent width'
		elif [ $p = vel ]
		then 
			cap='Velocity map'
		elif [ $p = sigma ]
		then
			cap='Velocity dispersion ($\mathrm{\sigma}$) map'
		elif [ $p = h3 ]
		then
			cap='Third Guass-Hermite moment (h3) map'
		elif [ $p = h4 ]
		then 
			cap='Fourth Guass-Hermite moment (h4) map'
		fi

		echo '            \caption{'${e^}$cap' for each galaxy in the sample. Plots are ordered roughly in peak stellar velocity}'
		echo '            \label{fig:'$c'_'$p'}'
		echo '        \end{figure*}'
		echo ''
		echo ''
	done
	echo ''
	echo ''
	echo ''
done

echo '\section{Kinematic Maps}'
echo '    \label{sec:kinmaps}'

for c in stellar gas
do
	if [ $c = stellar ]
	then 
		comp='vel sigma h3 h4'
		d=$c
	else
		comp='vel sigma'
		d='ionized gas'
	fi

	echo '    \subsection{'${d^}' kinematics maps}'
	echo '        \label{subsec:'$c'maps}'

	for p in $( echo $comp)
	do
		echo '        \begin{figure*}'
		echo '            \centering'
		# alphabetical 
		# for g in eso443-g024 ic1459 ic1531 ic4296 ngc0612 ngc1399 ngc3100 ngc3557 ngc7075 pks0718-34
		# rotational velocity
		for g in ngc0612 ngc3557 ngc3100 ic1459 pks0718-34 ic4296 ngc7075 ic1531 ngc1399 eso443-g024     
		do
			echo '            \includegraphics[width=0.245\textwidth]{plots/'$g'_'$c'_'$p'.png}'
		done

		if [ $p = vel ]
		then 
			cap='Velocity map'
		elif [ $p = sigma ]
		then
			cap='Velocity dispersion ($\mathrm{\sigma}$) map'
		elif [ $p = h3 ]
		then
			cap='Third Guass-Hermite moment (h3) map'
		elif [ $p = h4 ]
		then 
			cap='Fourth Guass-Hermite moment (h4) map'
		fi

		echo '            \caption{'${d^}' '$cap' for each galaxy in the sample. Plots are ordered roughly in peak stellar velocity}'
		echo '            \label{fig:'$c'_'$p'}'
		echo '        \end{figure*}'
		echo ''
		echo ''
	done
	echo ''
	echo ''
	echo ''
done