use std::collections::HashSet;
use std::path::Path;
use anyhow::Result;
use hpo2gene_parser::{GenePhenotypeLazyMapper, GenePhenotypeMapping}; 

#[test]
fn test_lazy_mapper_from_sample_file() -> Result<()> {
    let path = Path::new("sample_phenotype_to_genes.txt");
    let mapper = GenePhenotypeLazyMapper::from_file(path);

    let gene_to_hpo = mapper.gene_to_hpo()?;
    assert!(gene_to_hpo.contains_key("ALG14"));
    assert!(gene_to_hpo["ALG14"].contains("HP:0002188"));

    let hpo_to_genes = mapper.hpo_to_genes()?;
    assert!(hpo_to_genes["HP:0002188"].contains("ALG14"));
    assert!(hpo_to_genes["HP:0002188"].contains("AIFM1"));

    let gene_to_disease = mapper.gene_to_diseases()?;
    let diseases = &gene_to_disease["ALG14"];
    assert!(diseases.contains("OMIM:619031"));
    assert!(diseases.contains("OMIM:619036"));

    let hpo_name_to_ncbi = mapper.hpo_name_to_ncbi_ids()?;
    assert!(hpo_name_to_ncbi["Delayed CNS myelination"].contains("199857"));

    let ncbi_to_symbol = mapper.ncbi_to_symbol()?;
    assert_eq!(ncbi_to_symbol.get("199857"), Some(&"ALG14".to_string()));
    assert_eq!(ncbi_to_symbol.get("9131"), Some(&"AIFM1".to_string()));

    Ok(())
}

