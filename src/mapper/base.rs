use std::collections::{HashMap, HashSet};
use serde::Deserialize;
use serde::de::{self, Deserializer};
use anyhow::Result;
use ontolius::TermId;
use std::str::FromStr;

#[derive(Debug, Deserialize, Clone)]
pub struct PhenotypeRecord {
    #[serde(deserialize_with = "parse_term_id")]
    pub hpo_id: TermId,

    pub hpo_name: String,
    pub ncbi_gene_id: String,
    pub gene_symbol: String,

    #[serde(deserialize_with = "parse_term_id")]
    pub disease_id: TermId,
}


fn parse_term_id<'de, D>(deserializer: D) -> Result<TermId, D::Error>
where
    D: Deserializer<'de>,
{
    let s: &str = Deserialize::deserialize(deserializer)?;
    TermId::from_str(s).map_err(de::Error::custom)
}

#[derive(Debug, Clone, Copy)]
pub enum StringField {
    HpoName,
    NcbiGeneId,
    GeneSymbol,
}

#[derive(Debug, Clone, Copy)]
pub enum TermIdField {
    HpoId,
    DiseaseId,
}


pub trait GenePhenotypeMapping {
    /// Map gene symbol to HPO ID
    fn gene_to_hpo(&self) -> Result<HashMap<String, HashSet<TermId>>>;
    
    /// Map HPO ID to gene symbol
    fn hpo_to_genes(&self) -> Result<HashMap<TermId, HashSet<String>>>;

    /// Map gene symbol to disease IDs
    fn gene_to_diseases(&self) -> Result<HashMap<String, HashSet<TermId>>>;

    /// Map HPO name to NCBI gene ID
    fn hpo_name_to_ncbi_ids(&self) -> Result<HashMap<String, HashSet<String>>>;

    /// Map NCBI Gene ID to Gene Symbol
    fn ncbi_to_symbol(&self) ->  Result<HashMap<String, String>> ;
}
