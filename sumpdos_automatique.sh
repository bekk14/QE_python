#!/bin/bash

# ==============================================================================
#  Script pour extraire et sommer la PDOS pour chaque élément et chaque
#  orbitale atomique présente dans le répertoire de calcul.
# ==============================================================================

PREFIX=""
#upadate 
SUMPDOS="sumpdos.x" 

log() {
    echo "[INFO] $1"
}

# Détection automatique du préfixe si non défini
if [ -z "$PREFIX" ]; then
    first_file=$(find . -name "*.pdos_tot" -print -quit)
    if [ -n "$first_file" ]; then
        PREFIX=$(basename "$first_file" .pdos_tot)
        log "Préfixe détecté automatiquement : '${PREFIX}'"
    else
        echo "[ERREUR] Aucun fichier .pdos_tot trouvé. Impossible de deviner le préfixe."
        echo "Veuillez définir la variable PREFIX en haut du script."
        exit 1
    fi
fi
# Étape 1: Trouver tous les éléments uniques
log "Détection des éléments uniques..."
ELEMENTS=$(find . -name "${PREFIX}.pdos_atm#*" | sed -n 's/.*(\(.*\))_wfc.*/\1/p' | sort -u)
if [ -z "$ELEMENTS" ]; then
    echo "[ERREUR] Aucun élément détecté. Vérifiez que les fichiers PDOS existent."
    exit 1
fi
log "Éléments trouvés : $ELEMENTS"

# Étape 2: Trouver toutes les orbitales uniques
log "Détection des orbitales uniques..."
ORBITALS=$(find . -name "${PREFIX}.pdos_atm#*" | sed -n 's/.*_wfc#.*(\(.*\))/\1/p' | sort -u)
if [ -z "$ORBITALS" ]; then
    echo "[ERREUR] Aucune orbitale détectée. Vérifiez que les fichiers PDOS existent."
    exit 1
fi
log "Orbitales trouvées : $ORBITALS"

echo "----------------------------------------------------------------"
log "Démarrage du processus de sommation..."
echo "----------------------------------------------------------------"

# Étape 3: Boucler sur chaque élément et chaque orbitale pour exécuter sumpdos.x
for element in $ELEMENTS
do
    for orbital in $ORBITALS
    do
        search_pattern="${PREFIX}.pdos_atm#*\(${element}\)_wfc#*\(${orbital}\)"

        output_file="${element}_${orbital}_sum.dat"

        if compgen -G "$search_pattern" > /dev/null; then
            log "Traitement de : Élément=${element}, Orbitale=${orbital}"
            echo "  > Commande : sumpdos.x '$search_pattern' > $output_file"

            # Exécute la commande sumpdos.x
            $SUMPDOS $search_pattern > "$output_file"

            log "  > Fichier de sortie créé : ${output_file}"
        else
            log "Ignoré : Aucune orbitale '${orbital}' trouvée pour l'élément '${element}'."
        fi
    done
done

echo "----------------------------------------------------------------"
log "Processus de sommation terminé avec succès."
echo "----------------------------------------------------------------"

