# FERI K-MER Counter
Contador de k-mers en archivos .fasta, .fna o fastq. Básicamente es un contador que lee directamente las secuencias de DNA e ignora cualquier tipo de línea de calidad o comentario, o sea que trata a los archivos como simples repositorios de texto.

En este código se utiliza una tabla hash para guardar las secuencias, se ha probado de 10-mer hasta 30-mer en archivos no mayores a 800mb, obviamente el tiempo de lectura variará considerando los núcleos de procesamiento y condiciones del equipo que corra el código.
