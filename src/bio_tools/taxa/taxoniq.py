import taxoniq

def get_rank_name(taxon_or_name, rank_name="family"):
    """
    Return the scientific name of the requested rank for a taxon.

    Parameters
    ----------
    taxon_or_name : taxoniq.Taxon or str
        Either a Taxon object or a scientific name string (e.g. "Cicer arietinum").
    rank_name : str
        Rank to find (e.g. "family", "order", "genus").

    Returns
    -------
    str or None
        The scientific name at the requested rank, or None if not found.
    """
    # Accept either a Taxon object or a species name
    if isinstance(taxon_or_name, str):
        tax = taxoniq.Taxon(scientific_name=taxon_or_name)
    else:
        tax = taxon_or_name

    # 1) Preferred path: use ranked_lineage if available
    ranked = getattr(tax, "ranked_lineage", None)
    if ranked:
        for node in ranked:
            # many taxoniq versions store rank as an object with .name
            r = getattr(node, "rank", None)
            rname = None
            if r is not None:
                rname = getattr(r, "name", None) or getattr(r, "rank", None)
            else:
                # fallback: node.rank might be a plain string in some builds
                rname = getattr(node, "rank", None)
            if isinstance(rname, str) and rname.lower() == rank_name.lower():
                return node.scientific_name

    # 2) Fallback: walk up parent chain checking node.rank
    t = tax
    visited = set()
    while t is not None:
        # guard against cycles or repeated calls
        tid = getattr(t, "scientific_name", None) or id(t)
        if tid in visited:
            break
        visited.add(tid)

        r = getattr(t, "rank", None)
        rname = None
        if r is not None:
            rname = getattr(r, "name", None) or getattr(r, "rank", None)
        else:
            rname = getattr(t, "rank", None)

        if isinstance(rname, str) and rname.lower() == rank_name.lower():
            return getattr(t, "scientific_name", None)
        t = getattr(t, "parent", None)

    # not found
    return None
