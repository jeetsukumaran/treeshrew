
#include "reader.hpp"

namespace treeshrew {

TreeshrewBlockReader::TreeshrewBlockReader() :
        MultiFormatReader(-1, NxsReader::IGNORE_WARNINGS),
        reader_(nullptr) {
    this->SetWarningOutputLevel(AMBIGUOUS_CONTENT_WARNING);
    this->NCL_BLOCKTYPE_ATTR_NAME = "TREESHREW";
    NxsBlock::Reset();
    this->reader_ = this;
    NxsReader::Add(this);
    this->set_nxs_reader(this);
}

void TreeshrewBlockReader::read(const std::string& filepath, const std::string& format) {
    this->reader_->ReadFilepath(filepath.c_str(), MultiFormatReader::NEXUS_FORMAT);
    unsigned num_taxa_blocks = this->reader_->GetNumTaxaBlocks();
    std::cout << num_taxa_blocks << std::endl;
}

} // namespace treeshrew

