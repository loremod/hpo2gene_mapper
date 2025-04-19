mod base;
mod streaming_mapper;
mod memory_mapper;

pub use base::GenePhenotypeMapping;
pub use memory_mapper::GenePhenotypeMemoryMapper;
pub use streaming_mapper::GenePhenotypeStreamingMapper;
pub use base::PhenotypeRecord;