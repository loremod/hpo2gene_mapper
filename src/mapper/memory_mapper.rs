use std::collections::{HashMap, HashSet};
use std::fs::File;
use csv::ReaderBuilder;
use anyhow::Result;
use ontolius::TermId;
use std::path::Path;

use super::GenePhenotypeMapping;
use super::PhenotypeRecord;
use super::base::{StringField, TermIdField};

pub struct GenePhenotypeMemoryMapper {
    records: Vec<PhenotypeRecord>,
}

impl GenePhenotypeMapping for GenePhenotypeMemoryMapper{
    fn gene_to_hpo(&self) -> Result<HashMap<String, HashSet<TermId>>> {
        Ok(self.build_map_str_termid(StringField::GeneSymbol, TermIdField::HpoId))
    }

    fn hpo_to_genes(&self) -> Result<HashMap<TermId, HashSet<String>>> {
        Ok(self.build_map_termid_str(TermIdField::HpoId, StringField::GeneSymbol))
    }

    fn gene_to_diseases(&self) -> Result<HashMap<String, HashSet<TermId>>> {
        Ok(self.build_map_str_termid(StringField::GeneSymbol, TermIdField::DiseaseId))
    }

    fn hpo_name_to_ncbi_ids(&self) -> Result<HashMap<String, HashSet<String>>> {
        Ok(self.build_map_str_str(StringField::HpoName, StringField::NcbiGeneId))
    }

    fn ncbi_to_symbol(&self) ->  Result<HashMap<String, String>> {
        let mut map = HashMap::new();

        for record in &self.records {
            // Insert only if not already present
            map.entry(record.ncbi_gene_id.clone())
                .or_insert_with(|| record.gene_symbol.clone());
        }

       Ok(map)
    }
}



impl GenePhenotypeMemoryMapper {
    /// Constructor: Load dataset from file
    pub fn from_file<P: AsRef<Path>>(path: P) -> Result<Self> {
        let file = File::open(path)?;
        let mut rdr = ReaderBuilder::new().delimiter(b'\t').from_reader(file);

        let mut records = Vec::new();
        for result in rdr.deserialize() {
            let record: PhenotypeRecord = result?;
            records.push(record);
        }

        Ok(Self { records })
    }

    /// Internal method: Create a mapping from any string field to any string field
    fn build_map_str_str(&self, key: StringField, value: StringField) -> HashMap<String, HashSet<String>> {
        let mut map = HashMap::new();

        for record in &self.records {
            let k = Self::extract_string_field(record, key);
            let v = Self::extract_string_field(record, value);

            map.entry(k)
                .or_insert_with(HashSet::new)
                .insert(v);
        }

        map
    }

    /// Internal method: Create a mapping from any TermId to any string field
    fn build_map_termid_str(&self, key: TermIdField, value: StringField) -> HashMap<TermId, HashSet<String>> {
        let mut map = HashMap::new();

        for record in &self.records {
            let k = Self::extract_termid_field(record, key);
            let v = Self::extract_string_field(record, value);

            map.entry(k)
                .or_insert_with(HashSet::new)
                .insert(v);
        }

        map
    }

    /// Internal method: Create a mapping from any string to any TermId field
    fn build_map_str_termid(&self, key: StringField, value: TermIdField) -> HashMap<String, HashSet<TermId>> {
        let mut map = HashMap::new();

        for record in &self.records {
            let k = Self::extract_string_field(record, key);
            let v = Self::extract_termid_field(record, value);

            map.entry(k)
                .or_insert_with(HashSet::new)
                .insert(v);
        }

        map
    }

    /// Extract field value from a record, where the value is a String
    fn extract_string_field(record: &PhenotypeRecord, field: StringField) -> String {
        match field {
            StringField::HpoName => record.hpo_name.clone(),
            StringField::NcbiGeneId => record.ncbi_gene_id.clone(),
            StringField::GeneSymbol => record.gene_symbol.clone(),
        }
    }

    /// Extract field value from a record, where the value is a TermId
    fn extract_termid_field(record: &PhenotypeRecord, field: TermIdField) -> TermId {
        match field {
            TermIdField::HpoId => record.hpo_id.clone(),
            TermIdField::DiseaseId => record.disease_id.clone(),
        }
    }
}


#[cfg(test)]
mod tests {

    use super::*;

    #[test]
    fn test_memory_mapper_from_sample_file() -> Result<()> {
        let hpo1: TermId = "HP:0002188".parse().unwrap();
        
        let omim1: TermId = "OMIM:619031".parse().unwrap();
        let omim2: TermId = "OMIM:619036".parse().unwrap();


        let path = Path::new("data/sample_phenotype_to_genes.txt");
        let mapper = GenePhenotypeMemoryMapper::from_file(path).unwrap();

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