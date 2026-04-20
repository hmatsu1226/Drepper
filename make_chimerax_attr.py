#!/usr/bin/env python3
import sys
import csv
import matplotlib.pyplot as plt

def main():
    if len(sys.argv) != 5:
        print("Usage: python make_chimerax_attr.py <input.tsv> <protein_id> <output.defattr> <output.cxc>")
        sys.exit(1)

    input_tsv = sys.argv[1]
    target_id = sys.argv[2]
    output_attr = sys.argv[3]
    output_cxc = sys.argv[4]

    found = False

    with open(input_tsv, "r", encoding="utf-8") as f:
        reader = csv.reader(f, delimiter="\t")

        for row in reader:
            if not row:
                continue
            if row[0] != target_id:
                continue

            found = True

            cluster_str = row[-1].strip()
            cluster_ids = [int(x) for x in cluster_str.split(",") if x.strip() != ""]

            # --- attr file ---
            with open(output_attr, "w") as out:
                out.write("attribute: cluster\n")
                out.write("match mode: 1-to-1\n")
                out.write("recipient: residues\n\n")

                for i, cid in enumerate(cluster_ids, start=1):
                    out.write(f"\t:{i}\t{cid}\n")   # ← タブ重要

            print(f"[OK] wrote attr: {output_attr}")

            # --- color（tab20） ---
            unique_clusters = sorted(set(cluster_ids))
            cmap = plt.cm.get_cmap('tab20', len(unique_clusters))

            cluster_to_color = {}
            for i, cid in enumerate(unique_clusters):
                r, g, b, _ = cmap(i)
                r, g, b = int(r*255), int(g*255), int(b*255)
                cluster_to_color[cid] = (r, g, b)

            # --- cxc file ---
            with open(output_cxc, "w") as out:
                out.write("color #1F77B4\n")

                for cid, (r, g, b) in cluster_to_color.items():
                    hexcolor = f"#{r:02X}{g:02X}{b:02X}"
                    out.write(f"colordef c{cid} {hexcolor}\n")

                for cid in unique_clusters:
                    out.write(f"color ::cluster={cid} c{cid}\n")

            print(f"[OK] wrote cxc: {output_cxc}")

            return

    if not found:
        print(f"[ERROR] Protein ID not found: {target_id}")
        sys.exit(1)

if __name__ == "__main__":
    main()