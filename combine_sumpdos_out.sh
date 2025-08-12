#!/bin/bash

# ==============================================================================
#  Script pour combiner tous les fichiers PDOS sommés genere par sumpdos.x > (*_sum.dat) en un
#  seul fichier de données, avec les colonnes côte à côte.
#                 mis a jour 08/2025  
# !! apres sumpdos.x ..pdos_atm....  >  .._sum.dat 
# ==============================================================================

# --- Configuration ---
OUTPUT_FILE="pdos_total_combined.dat"

 
log() {
    echo "[INFO] $1"
}

# Vérifie si la commande 'paste' est disponible
if ! command -v paste &> /dev/null; then
    echo "[ERREUR] La commande 'paste' n'a pas été trouvée. Ce script ne peut pas continuer."
    exit 1
fi

# 
# On trie les fichiers pour avoir un ordre cohérent dans les colonnes
DATA_FILES=$(find . -maxdepth 1 -name "*_sum.dat" | sort)

if [ -z "$DATA_FILES" ]; then
    echo "[ERREUR] Aucun fichier de données '*_sum.dat' trouvé."
    echo "Veuillez d'abord lancer le script 'extract_all_pdos.sh'."
    exit 1
fi

log "Fichiers de données trouvés à combiner :"
echo "$DATA_FILES"
echo "----------------------------------------------------------------"

# --- Création de l'en-tête (header) pour le fichier de sortie ---
HEADER="Energy(eV)"

# Boucle sur chaque fichier pour construire le reste de l'en-tête
for file in $DATA_FILES; do
    # Extrait le nom de base (ex: Cs_s) du nom de fichier (ex: ./Cs_s_sum.dat)
    base_name=$(basename "$file" _sum.dat)
    HEADER="${HEADER}\t${base_name}"
done

# Écrit l'en-tête dans le fichier de sortie
# L'option -e permet d'interpréter le caractère de tabulation \t
echo -e "$HEADER" > "$OUTPUT_FILE"
log "En-tête créé pour le fichier de sortie : ${OUTPUT_FILE}"
log "En-tête : $HEADER"
echo "----------------------------------------------------------------"


# --- Combinaison des données ---


# Isole le premier fichier pour prendre sa colonne d'énergie
first_file=$(echo "$DATA_FILES" | head -n 1)


# Ex: <(cut -f2 -d' ' Cs_s_sum.dat) <(cut -f2 -d' ' K_p_sum.dat) ...
paste_commands=""
for file in $DATA_FILES; do
    # On utilise la substitution de processus <(...)
    # cut -f2 -d' ' : coupe la 2ème colonne, en utilisant l'espace comme délimiteur
    paste_commands+=" <(cut -f2 -d' ' \"$file\")"
done


log "Combinaison des données en cours... Veuillez patienter."
eval "paste <(cut -f1 -d' ' \"$first_file\") $paste_commands" >> "$OUTPUT_FILE"

log "Toutes les données ont été combinées avec succès dans le fichier : ${OUTPUT_FILE}"
echo "----------------------------------------------------------------"
echo "Vous pouvez maintenant utiliser ce fichier pour tracer vos graphiques."

