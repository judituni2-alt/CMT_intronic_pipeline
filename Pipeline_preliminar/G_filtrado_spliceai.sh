#!/bin/bash

set -euo pipefail


#------------------------------ARGUMENTOS---------------------------------


if [[ $# -lt 2 ]]; then
    echo "Uso: $0 input.vcf.gz salida_prefijo"
    exit 1
fi

VCF_IN="$1"
PREFIX="$2"
UMBRAL="${3:-0.5}"   # umbral opcional como 3er argumento; por defecto 0.5

TSV_RAW="${PREFIX}_spliceai.tsv"
TSV_FILTRADO="${PREFIX}_spliceai_${UMBRAL}.tsv"

if [[ ! -f "$VCF_IN" ]]; then
    echo "ERROR: no se encuentra el VCF de entrada: $VCF_IN"
    exit 1
fi

if ! command -v bcftools &> /dev/null; then
    echo "ERROR: 'bcftools' no está instalado o no está en el PATH."
    exit 1
fi


#-----------------EXTRAER CHROM/POS/REF/ALT/SpliceAI CON bcftools query----------------


echo "Extrayendo campos de $VCF_IN..."

# Añadimos una línea de cabecera manualmente, ya que 'bcftools query' no la
# genera por sí solo. Esto es lo que luego usará awk para saber qué columna
# es cuál (aunque en este caso, como SpliceAI viene en una sola columna con
# formato interno, el propio awk la vuelve a partir más abajo).
{
    printf 'CHROM\tPOS\tREF\tALT\tSpliceAI\n'
    bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/SpliceAI\n' "$VCF_IN"
} > "$TSV_RAW"

echo "Extracción completada: $TSV_RAW"

#---------------------FILTRAR POR DELTA SCORE >= UMBRAL-----------------

echo "Filtrando variantes con algún delta score >= $UMBRAL..."

awk -F'\t' -v umbral="$UMBRAL" '
BEGIN { OFS="\t" }
NR==1 { print; next }
{
    spliceai_field = $5
    incluir = 0

    # Puede haber varias anotaciones separadas por coma (multi-transcrito)
    n = split(spliceai_field, anotaciones, ",")

    for (i = 1; i <= n; i++) {
        # Cada anotación separada por "|":
        # 1=ALLELE 2=SYMBOL 3=DS_AG 4=DS_AL 5=DS_DG 6=DS_DL ...
        m = split(anotaciones[i], partes, "|")
        if (m < 6) continue   # anotación vacía o mal formada, se ignora

        ds_ag = partes[3] + 0
        ds_al = partes[4] + 0
        ds_dg = partes[5] + 0
        ds_dl = partes[6] + 0

        if (ds_ag >= umbral || ds_al >= umbral || ds_dg >= umbral || ds_dl >= umbral) {
            incluir = 1
            break
        }
    }

    if (incluir) print
}
' "$TSV_RAW" > "$TSV_FILTRADO"

N_TOTAL=$(($(wc -l < "$TSV_RAW") - 1))
N_FILTRADAS=$(($(wc -l < "$TSV_FILTRADO") - 1))

echo ""
echo "=========================================="
echo "Proceso completado."
echo "Variantes totales:        $N_TOTAL"
echo "Variantes con DS >= $UMBRAL: $N_FILTRADAS"
echo "Fichero completo:   $TSV_RAW"
echo "Fichero filtrado:   $TSV_FILTRADO"
echo "=========================================="
