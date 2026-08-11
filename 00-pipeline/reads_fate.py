#!/usr/bin/env python3

from pathlib import Path
import pandas as pd
import argparse
import re

PATTERN = re.compile(r'_[ATCG]+(?:_(?:negative|positive))?$')

def clean_id(id_string):
    return PATTERN.sub('', id_string)


def load_reads(files):
    """
    Charge les reads de chaque étape.
    """
    step_reads = {}

    for file in files:
        with open(file) as fh:
            reads = {
                clean_id(line.strip())
                for line in fh
                if line.strip()
            }

        step_reads[file.stem] = reads

    return step_reads


def build_presence_matrix(step_reads):
    """
    Construit une matrice présence/absence.
    """
    all_reads = set()

    for reads in step_reads.values():
        all_reads.update(reads)

    matrix = pd.DataFrame(index=sorted(all_reads))

    for step, reads in step_reads.items():
        matrix[step] = matrix.index.isin(reads).astype(int)

    matrix.index.name = "read_id"

    return matrix


def summarize_steps(step_reads):
    """
    Nombre de reads par étape.
    """
    summary = pd.DataFrame({
        "step": list(step_reads.keys()),
        "n_reads": [len(v) for v in step_reads.values()]
    })

    return summary



def main():

    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--input-dir",
        required=True,
        help="Répertoire contenant les fichiers txt"
    )

    parser.add_argument(
        "--output-dir",
        required=True,
        help="Répertoire de sortie"
    )

    args = parser.parse_args()

    input_dir = Path(args.input_dir)
    output_dir = Path(args.output_dir)

    output_dir.mkdir(parents=True, exist_ok=True)

    files = sorted(input_dir.glob("*.rnames"))

    if not files:
        raise ValueError("Aucun fichier TXT trouvé.")

    print(f"{len(files)} étapes détectées")

    step_reads = load_reads(files)

    # Matrice présence/absence
    presence = build_presence_matrix(step_reads)
    presence.to_csv(
        output_dir / "read_presence.tsv",
        sep="\t"
    )

    # Résumé
    summary = summarize_steps(step_reads)
    summary.to_csv(
        output_dir / "step_summary.tsv",
        sep="\t",
        index=False
    )

    print("Analyse terminée")


if __name__ == "__main__":
    main()