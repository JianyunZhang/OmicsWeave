"""
Multi-omics Cross-referencing Mappings
=======================================

This package provides two categories of mappings:

1. **Re-exported from ``omicsweave.core.idmap``** — canonical ID registries.
   New code should import directly from ``omicsweave.core.idmap``.

2. **Multi-omics specific mappings** — flux reaction → gene / protein / metabolite
   bridges loaded from YAML config files.  These are curated biological knowledge
   specific to the MFA-based multi-omics integration and are validated against the
   ``core.idmap`` registries at import time.
"""

# ===========================================================================
# 1. Re-exported from omicsweave.core.idmap
# ===========================================================================
from omicsweave.core.idmap.reaction_pathway_mapping import (
    REACTION_REGISTRY,
    get_reaction_info,
    get_kegg_reaction,
    get_kegg_pathway,
    get_reactome_id,
    get_ec_number,
    get_description,
    list_reactions_by_pathway,
    list_reactions_with_ec,
    get_all_pathways,
)
from omicsweave.core.idmap.mouse_protein_id_mapping import (
    MOUSE_PROTEIN_UNIPROT_MAPPING,
    get_uniprot_id,
    get_ec_for_protein,
)
from omicsweave.core.idmap.metabolite_id_mapping import (
    METABOLITE_DATABASE_IDS,
    lookup_metabolite_ids,
    annotate_metabolites,
    get_hmdb_id,
    get_kegg_compound,
    get_chebi_id,
)
from omicsweave.core.idmap.mouse_gene_id_mapping import (
    GENE_ANNOTATIONS,
    GENE_ALIAS_LOOKUP,
    resolve_gene_symbol,
    get_ncbi_gene_id,
    get_ensembl_gene_id,
    get_gene_name,
    get_gene_info,
)

# ===========================================================================
# 2. Multi-omics specific: flux → omics mappings  (YAML-backed)
# ===========================================================================
from .flux_gene_mapping import FLUX_GENE_MAPPING
from .flux_protein_mapping import (
    REACTION_ENZYME_MAPPING,
    get_all_related_proteins,
    get_related_proteins,
    get_protein_mapping_details,
)
from .flux_metabolite_mapping import (
    FLUX_METABOLITE_MAPPING,
    get_related_metabolites,
    get_metabolite_roles,
)
from .feature_resolver import (
    IdentifierMatch,
    MetaboliteRoleMatch,
    match_gene_features,
    match_metabolite_roles,
    match_protein_features,
)
from .audit import MappingIntegrityAudit, audit_mapping_integrity


# ===========================================================================
# 3. Convenience validation
# ===========================================================================

def validate_all_mappings() -> dict[str, list[str]]:
    """Run validation on all three flux mapping modules.

    Returns a dict mapping module name → list of warning strings.
    Modules that load cleanly produce an empty list.
    """
    from .flux_gene_mapping import _validate_genes
    from .flux_protein_mapping import _validate_proteins
    from .flux_metabolite_mapping import _validate_metabolites

    results = {
        "flux_gene_mapping": _validate_genes(FLUX_GENE_MAPPING),
        "flux_protein_mapping": _validate_proteins(REACTION_ENZYME_MAPPING),
        "flux_metabolite_mapping": _validate_metabolites(FLUX_METABOLITE_MAPPING),
    }

    for name, warns in results.items():
        if warns:
            print(f"[{name}] {len(warns)} warning(s):")
            for w in warns[:10]:
                print(f"  - {w}")
            if len(warns) > 10:
                print(f"  ... and {len(warns) - 10} more")
        else:
            print(f"[{name}] OK")

    return results


__all__ = [
    # --- Re-exported from common/ (prefer omicsweave.core) ---
    "REACTION_REGISTRY",
    "get_reaction_info",
    "get_kegg_reaction",
    "get_kegg_pathway",
    "get_reactome_id",
    "get_ec_number",
    "get_description",
    "list_reactions_by_pathway",
    "list_reactions_with_ec",
    "get_all_pathways",
    "MOUSE_PROTEIN_UNIPROT_MAPPING",
    "get_uniprot_id",
    "get_ec_for_protein",
    "METABOLITE_DATABASE_IDS",
    "lookup_metabolite_ids",
    "annotate_metabolites",
    "get_hmdb_id",
    "get_kegg_compound",
    "get_chebi_id",
    "GENE_ANNOTATIONS",
    "GENE_ALIAS_LOOKUP",
    "resolve_gene_symbol",
    "get_ncbi_gene_id",
    "get_ensembl_gene_id",
    "get_gene_name",
    "get_gene_info",
    # --- Multi-omics specific (YAML-backed) ---
    "FLUX_GENE_MAPPING",
    "REACTION_ENZYME_MAPPING",
    "FLUX_METABOLITE_MAPPING",
    "get_all_related_proteins",
    "get_related_proteins",
    "get_protein_mapping_details",
    "get_related_metabolites",
    "get_metabolite_roles",
    "IdentifierMatch",
    "MetaboliteRoleMatch",
    "match_gene_features",
    "match_metabolite_roles",
    "match_protein_features",
    "MappingIntegrityAudit",
    "audit_mapping_integrity",
    # --- Convenience ---
    "validate_all_mappings",
]
