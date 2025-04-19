use std::collections::{HashMap, HashSet};
use std::fs::File;
use std::path::{Path, PathBuf};
use csv::ReaderBuilder;
use anyhow::Result;
use ontolius::TermId;

// use ontolius::TermId;

use super::GenePhenotypeMapping;
use super::PhenotypeRecord;

pub struct GenePhenotypeStreamingMapper {
    path: PathBuf,
}

impl GenePhenotypeStreamingMapper {
    pub fn from_file<P: AsRef<Path>>(path: P) -> Result<Self> {
        Ok(Self {
            path: path.as_ref().to_path_buf(),
        })
    }

    fn parse_file(&self) -> Result<impl Iterator<Item = PhenotypeRecord>> {
        let file = File::open(&self.path)?;
        let rdr = ReaderBuilder::new().delimiter(b'\t').from_reader(file);
        Ok(rdr.into_deserialize::<PhenotypeRecord>().filter_map(Result::ok))
    }
}

impl GenePhenotypeMapping for GenePhenotypeStreamingMapper {
    fn gene_to_hpo(&self) -> Result<HashMap<String, HashSet<TermId>>> {
        let mut map: HashMap<String, HashSet<TermId>> = HashMap::new();
        for record in self.parse_file()? {
            map.entry(record.gene_symbol).or_default().insert(record.hpo_id);
        }
        Ok(map)
    }

    fn hpo_to_genes(&self) -> Result<HashMap<TermId, HashSet<String>>> {
        let mut map: HashMap<TermId, HashSet<String>> = HashMap::new();
        for record in self.parse_file()? {
            map.entry(record.hpo_id).or_default().insert(record.gene_symbol);
        }
        Ok(map)
    }

    fn gene_to_diseases(&self) -> Result<HashMap<String, HashSet<TermId>>> {
        let mut map: HashMap<String, HashSet<TermId>> = HashMap::new();
        for record in self.parse_file()? {
            map.entry(record.gene_symbol).or_default().insert(record.disease_id);
        }
        Ok(map)
    }

    fn hpo_name_to_ncbi_ids(&self) -> Result<HashMap<String, HashSet<String>>> {
        let mut map: HashMap<String, HashSet<String>> = HashMap::new();
        for record in self.parse_file()? {
            map.entry(record.hpo_name).or_default().insert(record.ncbi_gene_id);
        }
        Ok(map)
    }

    fn ncbi_to_symbol(&self) -> Result<HashMap<String, String>> {
        let mut map: HashMap<String, String> = HashMap::new();
        for record in self.parse_file()? {
            map.entry(record.ncbi_gene_id.clone())
                .or_insert(record.gene_symbol.clone());
        }
        Ok(map)
    }
}

#[cfg(test)]
mod tests {

    use super::*;


#[test]
fn test_streaming_mapper_from_sample_file() -> Result<()> {
    let hpo1: TermId = "HP:0002188".parse().unwrap();
    
    let omim1: TermId = "OMIM:619031".parse().unwrap();
    let omim2: TermId = "OMIM:619036".parse().unwrap();


    let path = Path::new("data/sample_phenotype_to_genes.txt");
    let mapper = GenePhenotypeStreamingMapper::from_file(path).unwrap();

    let gene_to_hpo = mapper.gene_to_hpo()?;
    assert!(gene_to_hpo.contains_key("ALG14"));
    assert!(gene_to_hpo["ALG14"].contains(&hpo1));

    let hpo_to_genes = mapper.hpo_to_genes()?;
    assert!(hpo_to_genes[&hpo1].contains("ALG14"));
    assert!(hpo_to_genes[&hpo1].contains("AIFM1"));

    let gene_to_disease = mapper.gene_to_diseases()?;
    let diseases = &gene_to_disease["ALG14"];
    assert!(diseases.contains(&omim1));
    assert!(diseases.contains(&omim2));

    let hpo_name_to_ncbi = mapper.hpo_name_to_ncbi_ids()?;
    assert!(hpo_name_to_ncbi["Delayed CNS myelination"].contains("199857"));

    let ncbi_to_symbol = mapper.ncbi_to_symbol()?;
    assert_eq!(ncbi_to_symbol.get("199857"), Some(&"ALG14".to_string()));
    assert_eq!(ncbi_to_symbol.get("9131"), Some(&"AIFM1".to_string()));

    Ok(())
}

}