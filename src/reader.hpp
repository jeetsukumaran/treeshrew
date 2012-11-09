#ifndef TREESHREW_READER_HPP
#define TREESHREW_READER_HPP

#include "ncl/nxsmultiformat.h"

namespace treeshrew {

class TreeshrewBlockReader : public NxsBlock, public MultiFormatReader {

    public:
        TreeshrewBlockReader();
        void read(const std::string& filepath, const std::string& format);
        void set_nxs_reader(MultiFormatReader * r) {
            this->reader_ = r;
        }

    private:
        MultiFormatReader *     reader_;


}; // TreeShrewBlockReader

} // namespace treeshrew

#endif
