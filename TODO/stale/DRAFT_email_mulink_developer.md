# Draft email — mulink developer

To: Lucas Diedrich `<diedrich@biochem.mpg.de>`  
Subject: APB and mulink — representing proteomics hierarchies in MuData

Dear Lucas,

I am developing the **anndata-proteomics-bridge (APB)**:
https://github.com/anndata-omics-bridge/anndata-proteomics-bridge

APB converts proteomics software output into AnnData and MuData while preserving the original
measurements, metadata, and provenance. Depending on what a vendor report provides, APB represents
fragment, ion/precursor, peptidoform, peptide, and protein quantification as separate modalities.

We recently reviewed mulink and the work initiated at the scverse proteomics workshop. It addresses
exactly the question we now face: how to represent the relationships between proteomics feature
levels in a scverse-native way. We would like APB to align with and use **mulink**, rather than
introduce a separate APB-specific linking model.

I would be very interested in your view on the intended conventions for proteomics hierarchies,
particularly:

- how fragment, precursor, peptide, protein, and gene features should be linked;
- the intended direction and meaning of those relationships;
- how mapping provenance should be represented, for example vendor-reported protein assignments
  versus peptide-to-protein matches derived from FASTA;
- how APB could produce mulink-compatible MuData without duplicating functionality;
- whether a small shared APB/mulink example would be useful.

Would you be open to discussing this by email or in a short call? I would be happy to share a small
public APB-generated MuData example as a concrete starting point.

Best regards,

Witold Wolski
