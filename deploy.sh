#!/usr/bin/env bash
set -euo pipefail

DISTRICTS=(queenstown bishan outram bukit-merah newton)

echo "=== Regenerating processed layers ==="

if [ -d data/data-gov-sg ]; then
  for d in "${DISTRICTS[@]}"; do
    echo ""
    echo "--- $d ---"
    python3 src/generate_subzone_summary.py --district "$d"
    python3 src/filter_district_layers.py --district "$d"
  done
else
  echo "WARNING: data/data-gov-sg/ not found, skipping data regeneration"
  echo "  Committed files in docs/geo/ will be used as-is"
fi

echo ""
echo "=== Pushing to main ==="
git push origin main
echo "=== Done ==="
