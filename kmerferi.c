// ---------------------------------------------------------------------------------------
// Autor: Felipe Alexander Correa Rodríguez
// Fecha: 2023-10-02
// Descripción: Este programa calcula la frecuencia de K-mers en un archivo FASTA y FASTQ.
// ---------------------------------------------------------------------------------------

// Librerías
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <omp.h>
#include <stdint.h>
#include <time.h>

// Definiciones de constantes
#define LONGITUD_KMERO 10  // Longitud del K-mero
#define MAX_SEQ 100000 // Número máximo de secuencias de K-meros, 100000 K-meros
int contador = 0; // contador de secuencias de ADN

// Esta estructura representa una secuencia de ADN y su frecuencia
typedef struct 
{
    char seq[LONGITUD_KMERO + 1]; 
    int freq;
} dnaSeqFreq;

typedef struct dnaSeqFreqNode 
{
    dnaSeqFreq *data;
    struct dnaSeqFreqNode *next;
} dnaSeqFreqNode;

dnaSeqFreqNode* hashtable[MAX_SEQ];

// Función hash, básicamente una función de dispersión de cadenas
uint32_t hash(char *str)  // variables de entrada: str, variable de salida: hash
{
    uint32_t hash = 0;
    int c;
    while ((c = *str++))
        hash = c + (hash << 6) + (hash << 16) - hash;
    return hash % MAX_SEQ;
}

// Función de búsqueda
dnaSeqFreq* buscar(char *seq) // variables de entrada: seq, variable de salida: node
{
    dnaSeqFreqNode* node = hashtable[hash(seq)];
    while (node != NULL && strcmp(node->data->seq, seq))
        node = node->next;
    return node ? node->data : NULL;
}

// Función de inserción, inserta una secuencia de ADN en la tabla hash
void insertar(char *seq)  // variables de entrada: seq, variable de salida: node
{
    dnaSeqFreq *s = buscar(seq);  // punto de entrada a la función buscar
    if (s == NULL) {
        uint32_t h = hash(seq);
        dnaSeqFreqNode *node = (dnaSeqFreqNode*)malloc(sizeof(dnaSeqFreqNode)); // asigna memoria
        dnaSeqFreq *s = (dnaSeqFreq*)malloc(sizeof(dnaSeqFreq)); // asigna memoria
        strncpy(s->seq, seq, LONGITUD_KMERO); // copia la secuencia de ADN
        s->seq[LONGITUD_KMERO] = '\0';  // termina la secuencia de ADN
        s->freq = 1; // inicializa la frecuencia
        node->data = s; // asigna la secuencia de ADN
        node->next = hashtable[h];
        hashtable[h] = node; // asigna el nodo
        contador++; // Incrementa el contador de secuencias de ADN
    } else {
        s->freq++;
    }
}

// Función principal
int main() {
    double tiempo_inicio = omp_get_wtime(); // tiempo de inicio
    int num_nucleos = omp_get_num_procs(); // número de nucleos
    FILE* archivo = fopen("ecoli.fasta", "r"); // abre el archivo de entrada

    if (archivo == NULL) 
    {
        printf("No se pudo abrir el archivo.\n");
        return 1;
    }

    char linea[256];

    while (fgets(linea, sizeof(linea), archivo)) 
    {
        // Ignoramos las líneas que no son secuencias. Esto para archivos FASTA y FASTQ.
        if (linea[0] == '@' || linea[0] == '>' || linea[0] == '+' || linea[0] == '#')
            continue; // continua con la siguiente iteración del ciclo while
        
        int stop = strlen(linea) - LONGITUD_KMERO + 1;
#pragma omp parallel for // paraleliza el ciclo for
        for (int i = 0; i < stop; i++) 
        {
            char seq[LONGITUD_KMERO+1] = {0};
            strncpy(seq, &linea[i], LONGITUD_KMERO);
            #pragma omp critical // sección crítica, solo un hilo a la vez
            {
                insertar(seq);
            }
        }
    }
    fclose(archivo);

    // Contar frecuencias de K-mers
    int freq_kmer[MAX_SEQ] = {0};
    for(int i = 0; i < MAX_SEQ; i++)
    {
        dnaSeqFreqNode* node = hashtable[i];
        while(node) 
        {
            dnaSeqFreq* s = node->data;
            if(s) 
            {
                freq_kmer[s->freq < 100 ? s->freq : 100]++;
            }
            node = node->next;
        }
    }
    
    FILE* archivo_salida = fopen("ecoli3.out", "w"); // Abre el archivo de salida
    // En este archivo se guardara la frecuencia de cada K-mer donde la primera columna es el K-mer y la segunda columna es la frecuencia
    if (archivo_salida == NULL) 
    {
        printf("No se pudo abrir el archivo de salida.\n");
        return 1;
    }

    for(int i = 1; i < MAX_SEQ; i++) // Recorre la tabla de hash
    {
        if(freq_kmer[i] > 0)
        {
            fprintf(archivo_salida, "%d -> %d\n", i, freq_kmer[i]);
        }
    }

    fclose(archivo_salida);

    double tiempo_fin = omp_get_wtime(); // pin de tiempo

// Impresión de resultados:
    printf("Calculo de secuencias de ADN en %d-mer finalizado.\n", LONGITUD_KMERO);
    printf("Numero de K-mers calculados: %d\n", contador);
    printf("Tiempo de ejecucion: %f segundos.\n", tiempo_fin - tiempo_inicio);
    printf("Numero de nucleos utilizados: %d\n", num_nucleos);

    return 0;
}


// ---------------------------------------------------------------------------------------
// Documentacion basica:
// Este codigo realiza el siguiente camino para contar k-mers:
// 1. Abre el archivo de entrada.
// 2. Lee el archivo de entrada.
// 3. Ignora las lineas que no son secuencias.
// 4. Crea un ciclo for para recorrer las secuencias.
// 5. Crea una sección crítica para insertar las secuencias en la tabla hash.
// (Una tabla hash es una estructura de datos que permite almacenar y recuperar datos,
// se ocupa aquí para contar las secuencias de ADN, en este caso, los k-mers)
// 6. Cierra el archivo de entrada.
// 7. Abre el archivo de salida.
// 8. Recorre la tabla hash.
// 9. Cierra el archivo de salida.
// 10. Imprime los resultados.
// ---------------------------------------------------------------------------------------
// Compilar y ejecutar:
// gcc -Wall -fopenmp -o kmerferi kmerferi.c
// kmerferi.exe