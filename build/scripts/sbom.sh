#
# SPDX-License-Identifier: MIT
# Copyright (c) "2025" . The DeepCausality Authors and Contributors. All Rights Reserved.
#

set -o errexit
set -o nounset
set -o pipefail

CRATES=(
    "data_manager"
    "deep_causality_physics"
    "experiments"
)

for CRATE_NAME in "${CRATES[@]}"; do
    echo "Generating SBOM for crate: $CRATE_NAME"

    if ! cargo sbom --cargo-package "$CRATE_NAME" --output-format=spdx_json_2_3 > "$CRATE_NAME"/"$CRATE_NAME"_sbom.spdx.json
    then
        echo "Failed to generate SBOM for $CRATE_NAME"
    fi

     if ! sha256sum "$CRATE_NAME"/"$CRATE_NAME"_sbom.spdx.json > "$CRATE_NAME"/"$CRATE_NAME"_sbom.spdx.json.sha
     then
        echo "Failed to generate HASH over SBOM for $CRATE_NAME"
     fi

done

echo "SBOM generation complete."
